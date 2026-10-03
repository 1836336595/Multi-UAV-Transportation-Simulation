function dump_sim_data(outDir)
%DUMP_SIM_DATA  跑一次仿真，把关键日志导出成 CSV + 诊断报告，供离线排查。
%
% 用法（工作目录必须在 git 目录下，或已 addpath 该目录）：
%     dump_sim_data();                       % 输出到 ..\_debug_out
%     dump_sim_data('F:\tmp\dump');          % 指定输出目录
%
% 产出（outDir 下）：
%     dump_report.txt       诊断报告（人读）：配置摘要 / 状态机时序 / 交接窗口峰值 /
%                           抬升段过冲与稳态 / 全程极值 / summary 全字段
%     dump_trajectory.csv   time, 模式码, tensionScale, 实际位置, 期望位置, 高度, 误差
%     dump_velocity.csv     负载速度
%     dump_tension.csv      实际张力, 期望张力, 推力占比
%     dump_link.csv         绳向误差(deg), 绳长, 松弛量, 绳向竖直流向分量
%     dump_attitude.csv     负载姿态误差, 负载角速度, 期望/实际偏航
%     dump_vehicle.csv      三机位置
%
% ★★ 本脚本是**调试工具**，不参与仿真，不改动任何仿真代码，只读 sim 的日志。
%    提交 GitHub 前可直接删除（或写进 .gitignore）。
% ★★ 关键设计：**期望位置不靠重新调用参考函数**（匿名句柄多输出调用在部分
%    MATLAB 版本会失败，且 TAUT_RAMP 期间的期望被参考平滑改写过，再算一遍会失真）。
%    这里用日志里的恒等式精确反解：
%        positionErrorVectorLog = 负载实际位置 - 控制器看到的期望位置
%    （见 crazyflie_slung_controller.m 中 ex = state.loadPosition - desired.position，
%      以及 simulation.m 里 sim.positionErrorVectorLog(:, k) = command.positionError）
%
% ★ 绳向误差的物理含义：controller.m 里 eqi = cross(q_id, q_i)，
%   所以 sim.linkErrorLog = ||cross|| = sin(夹角)。本脚本统一换算成**度**。

if nargin < 1 || isempty(outDir)
    thisDir = fileparts(mfilename('fullpath'));
    outDir = fullfile(thisDir, '..', '_debug_out');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outDir = char(outDir);
fprintf('[dump] 输出目录：%s\n', outDir);

% ------------------------------------------------------------------ 配置
cfg = crazyflie_slung_parameters();
if isfield(cfg, 'visualization')
    if isfield(cfg.visualization, 'plot'),    cfg.visualization.plot = false;    end
    if isfield(cfg.visualization, 'animate'), cfg.visualization.animate = false; end
end
% ★★ persistent 状态在同一个 MATLAB 会话里**不会**自动重置：如果本脚本连跑两次，
%    第二次会沿用上一次的 firedKinds / 状态，症状是"这次没发散"的假象
%    （本项目为此误判过一次）。必须显式 clear。
clear('crazyflie_slung_dynamics');

fprintf('[dump] 开始仿真 ...\n');
wallClock = tic;
sim = crazyflie_slung_simulation(cfg);
fprintf('[dump] 仿真结束：耗时 %.2f s，步数 %d\n', toc(wallClock), numel(sim.time));

% ------------------------------------------------------------------ 取日志
t = sim.time(:).';
nSteps = numel(t);
n = cfg.vehicle.count;

loadPos  = sim.loadPositionLog;             % 3 x nSteps
loadVel  = sim.loadVelocityLog;             % 3 x nSteps
errVec   = sim.positionErrorVectorLog;      % 3 x nSteps（实际 - 期望）
desPos   = loadPos - errVec;                % 反解出的期望位置（控制器真正看到的）
height   = -loadPos(3, :);                  % z 向下为正 ⇒ 高度取负
desH     = -desPos(3, :);
linkErr  = sim.linkErrorLog;                % n x nSteps，= sin(绳向夹角)
linkErrDeg = asind(min(max(linkErr, -1), 1));
T        = sim.tensionLog;                  % n x nSteps
% 期望张力：desiredTensionLog 是 3 x n x nSteps，取每个绳的模长 -> n x nSteps
dT       = sim.desiredTensionLog;           % 3 x n x nSteps
desT     = reshape(sqrt(sum(dT .^ 2, 1)), n, nSteps);
thrustPct = sim.thrustPctLog;               % n x nSteps
% 绳向单位向量（3 x n x nSteps）→ 竖直流向分量（悬停时 ≈ +1）
qNorm    = reshape(sqrt(sum(sim.linkUnitLog .^ 2, 1)), n, nSteps);
qVert    = reshape(sim.linkUnitLog(3, :, :), n, nSteps);
loadRate = sim.loadBodyRateLog;             % 3 x nSteps
yawRefDeg = rad2deg(sim.loadYawRefLog);
yawDeg    = rad2deg(sim.loadYawLog);
yawErrDeg = rad2deg(sim.loadYawErrorLog);
modeSeq  = sim.takeoffModeLog;              % 1 x nSteps
tScale   = sim.tensionScaleLog;             % 1 x nSteps

% 各绳最大绳向误差、各机最大推力占比（1 x nSteps）
linkErrWorst = max(linkErrDeg, [], 1);
thrustWorst  = max(thrustPct, [], 1);
errNorm      = sqrt(sum(errVec .^ 2, 1));

% 状态机时序（modeCode：0=SLACK 1=TAKEUP 2=TAUT_RAMP 3=ACTIVE
%                        4=LANDING_TAUT 5=LANDING_RELEASE）
modeNames = {'SLACK', 'TAKEUP', 'TAUT_RAMP', 'ACTIVE', ...
             'LANDING_TAUT', 'LANDING_RELEASE'};
enterTime = nan(1, 6);
leaveTime = nan(1, 6);
for c = 0:5
    idxFirst = find(modeSeq == c, 1, 'first');
    if ~isempty(idxFirst), enterTime(c + 1) = t(idxFirst); end
    idxLast = find(modeSeq == c, 1, 'last');
    if ~isempty(idxLast), leaveTime(c + 1) = t(idxLast); end
end

tRamp = enterTime(3);                        % 进入 TAUT_RAMP = 交接时刻

% 稳态窗口：与 simulation.m 的 summary 保持一致（max(1,round(0.8*nSteps)):nSteps）
steadyIdx = max(1, round(0.8 * nSteps)):nSteps;
% 只有绷紧段（模式 2/3/4）的稳态窗口才代表"悬挂运输"
tautMask = (modeSeq == 2) | (modeSeq == 3) | (modeSeq == 4);
steadyTautIdx = steadyIdx(tautMask(steadyIdx));
if isempty(steadyTautIdx)
    steadyTautIdx = find(tautMask);
end
% ★★ ACTIVE 段才是"悬挂运输"真正的稳态窗口。
%   末 20% 窗口会**包含降落段**（负载已经在下降/落地），报出来的"稳态误差"
%   严重误导。实测（2026-09-28 转储）：末 20% 给"高度 0.2322 m、位置误差 0.1435 m"，
%   看起来像存在稳态偏差；而 ACTIVE 段实际是"高度 0.35005 m（偏差 0.05 mm）、
%   位置误差均值 0.11 mm"。⇒ 稳态指标一律以 activeIdx 为准。
activeIdx = find(modeSeq == 3);
if isempty(activeIdx)
    activeIdx = steadyTautIdx;
end

targetH = -cfg.target.position(3);           % 目标离地高度 [m]

% ------------------------------------------------------------------ 报告
R = {};
R = pushLine(R, '==============================================================================');
R = pushLine(R, ' 仿真数据转储报告  —  %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
R = pushLine(R, '==============================================================================');

R = pushLine(R, '');
R = pushLine(R, '【A】配置摘要');
R = pushLine(R, '  四旋翼数量 n            : %d', n);
R = pushLine(R, '  负载质量 m0             : %.6g kg', cfg.payload.mass);
R = pushLine(R, '  负载尺寸 (a,b,c)        : [%s] m', numVecText(cfg.payload.size(:).', '%.4g '));
R = pushLine(R, '  单机质量 / 绳长         : %.6g kg / %.6g m', cfg.vehicle.mass, cfg.link.length);
R = pushLine(R, '  仿真时长 / dt / 步数    : %.4g s / %.6g s / %d', cfg.simulation.duration, cfg.simulation.dt, nSteps);
R = pushLine(R, '  参考轨迹                : 静态目标（定高工况，无外部轨迹函数）');
R = pushLine(R, '  target.position         : [%s] m  ⇒ 目标高度 %.6g m', ...
    numVecText(cfg.target.position(:).', '%.6g '), targetH);
if isfield(cfg, 'takeoff') && isfield(cfg.takeoff, 'enabled')
    R = pushLine(R, '  takeoff.enabled         : %d   landingEnabled = %d', ...
        cfg.takeoff.enabled, getOr(cfg.takeoff, 'landingEnabled', NaN));
    R = pushLine(R, '  independentHoverHeight   : %.6g m', getOr(cfg.takeoff, 'independentHoverHeight', NaN));
    R = pushLine(R, '  takeupDuration          : %.6g s', getOr(cfg.takeoff, 'takeupDuration', NaN));
    R = pushLine(R, '  preTensionSlack         : %.6g m', getOr(cfg.takeoff, 'preTensionSlack', NaN));
    R = pushLine(R, '  epsilonOn / epsilonOff  : %.6g / %.6g m', ...
        getOr(cfg.takeoff, 'epsilonOn', NaN), getOr(cfg.takeoff, 'epsilonOff', NaN));
    R = pushLine(R, '  confirmTime             : %.6g s', getOr(cfg.takeoff, 'confirmTime', NaN));
    R = pushLine(R, '  tensionRampTime         : %.6g s', getOr(cfg.takeoff, 'tensionRampTime', NaN));
    R = pushLine(R, '  landingDuration         : %.6g s', getOr(cfg.takeoff, 'landingDuration', NaN));
end
R = pushLine(R, '  loadController.kR       : [%s]', numVecText(cfg.loadController.kR(:).', '%.6g '));
R = pushLine(R, '  loadController.kOmega   : [%s]', numVecText(cfg.loadController.kOmega(:).', '%.6g '));
R = pushLine(R, '  kx / kv / ki            : [%s] / [%s] / [%s]', ...
    numVecText(cfg.loadController.kx(:).', '%.4g '), ...
    numVecText(cfg.loadController.kv(:).', '%.4g '), ...
    numVecText(cfg.loadController.ki(:).', '%.4g '));

R = pushLine(R, '');
R = pushLine(R, '【B】状态机时序（modeCode 0..5）');
for c = 0:5
    if ~isnan(enterTime(c + 1))
        R = pushLine(R, '  [%d] %-18s 进入 %8.4f s   离开 %8.4f s   持续 %8.4f s', ...
            c, modeNames{c + 1}, enterTime(c + 1), leaveTime(c + 1), ...
            leaveTime(c + 1) - enterTime(c + 1));
    end
end
takeoffOn = isfield(cfg, 'takeoff') && isfield(cfg.takeoff, 'enabled') ...
    && logical(cfg.takeoff.enabled);
if ~takeoffOn
    R = pushLine(R, '  ★ takeoff.enabled = false（或缺少该字段）：没有起飞状态机，');
    R = pushLine(R, '    全程直接走绷紧动力学，"交接"指标（C/D 两节）仅供参照。');
elseif isnan(tRamp)
    R = pushLine(R, '  ★ 从未进入 TAUT_RAMP —— 绳没有绷紧，后面所有"交接"指标无意义。');
else
    R = pushLine(R, '  ★ 交接时刻（进入 TAUT_RAMP）t = %.4f s', tRamp);
end

R = pushLine(R, '');
R = pushLine(R, '【C】交接窗口（[tRamp-0.5, tRamp+6] s）峰值');
if isnan(tRamp)
    w0 = t(1); w1 = t(end);
else
    w0 = max(t(1), tRamp - 0.5);
    w1 = min(t(end), tRamp + 6);
end
R = pushLine(R, '  窗口                   : [%.4f, %.4f] s', w0, w1);
[v, tv] = peakInWindow(linkErrWorst, t, w0, w1);
R = pushLine(R, '  绳向误差峰值（取最差绳）: %8.4f deg  @ t = %.4f s', v, tv);
[v, tv] = peakInWindow(height, t, w0, w1);
R = pushLine(R, '  负载高度峰值            : %8.4f m    @ t = %.4f s', v, tv);
[v, tv] = peakInWindow(errNorm, t, w0, w1);
R = pushLine(R, '  位置误差峰值            : %8.4f m    @ t = %.4f s', v, tv);
[v, tv, iL] = peakInWindowMatrix(T, t, w0, w1, 'min');
R = pushLine(R, '  最小张力                : %8.4f N    @ t = %.4f s（第 %d 根绳）', v, tv, iL);
[v, tv] = peakInWindow(thrustWorst, t, w0, w1);
R = pushLine(R, '  最大推力占比            : %8.4f %%   @ t = %.4f s', v, tv);
[v, tv] = peakInWindow(sqrt(sum(loadRate .^ 2, 1)), t, w0, w1);
R = pushLine(R, '  负载角速度峰值          : %8.4f rad/s @ t = %.4f s', v, tv);
[v, tv] = peakInWindow(abs(yawErrDeg), t, w0, w1);
R = pushLine(R, '  偏航跟踪误差峰值        : %8.4f deg  @ t = %.4f s', v, tv);
[v, tv] = peakInWindow(abs(tScale), t, w0, w1);
R = pushLine(R, '  tensionScale 峰值       : %8.4f      @ t = %.4f s', v, tv);

R = pushLine(R, '');
R = pushLine(R, '【D】抬升段：过冲 / 稳态偏差（直接回答"是过冲还是稳态偏差"）');
R = pushLine(R, '  目标高度                : %8.4f m', targetH);
if isnan(tRamp)
    liftW0 = t(1);
else
    liftW0 = tRamp;
end
liftIdx = find(t >= liftW0);
if ~isempty(liftIdx)
    [hMax, j] = max(height(liftIdx));
    R = pushLine(R, '  抬升段高度峰值          : %8.4f m    @ t = %.4f s', hMax, t(liftIdx(j)));
    R = pushLine(R, '  ⇒ 超出目标              : %+8.4f m', hMax - targetH);
else
    R = pushLine(R, '  抬升段：无采样');
end
if ~isnan(tRamp)
    for dtProbe = [1, 2, 5, 10, 20]
        tp = tRamp + dtProbe;
        if tp <= t(end)
            hp = interp1(t, height, tp, 'linear');
            R = pushLine(R, '  高度 @ tRamp%+3d s      : %8.4f m（偏差 %+8.4f m）', ...
                dtProbe, hp, hp - targetH);
        end
    end
end
R = pushLine(R, '  稳态高度（ACTIVE 段均值） : %8.4f m（偏差 %+8.4f m）← 真稳态', ...
    mean(height(activeIdx)), mean(height(activeIdx)) - targetH);
R = pushLine(R, '  稳态高度（末 20%% 绷紧段）: %8.4f m（偏差 %+8.4f m）', ...
    mean(height(steadyTautIdx)), mean(height(steadyTautIdx)) - targetH);
R = pushLine(R, '  稳态高度（末 20%% 全部）  : %8.4f m（偏差 %+8.4f m）← 含降落段，勿用', ...
    mean(height(steadyIdx)), mean(height(steadyIdx)) - targetH);
R = pushLine(R, '  ★ 判据：抬升段峰值偏差大、而 ACTIVE 段偏差 ≈ 0 ⇒ 是**过冲/瞬态**；');
R = pushLine(R, '    若 ACTIVE 段偏差也大 ⇒ 才是**稳态偏差**，该看积分器/重力前馈而不是阻尼。');
R = pushLine(R, '  ★ 末 20%% 窗口含降落段（负载已落地）⇒ 它的"稳态"数字必然偏大，只作参考。');

R = pushLine(R, '');
R = pushLine(R, '【E】全程极值');
[v, j] = max(errNorm);
R = pushLine(R, '  位置误差最大            : %8.4f m    @ t = %.4f s   分量 [%s]', ...
    v, t(j), numVecText(errVec(:, j).', '%+.4f '));
[v, j] = max(linkErrWorst);
R = pushLine(R, '  绳向误差最大（最差绳）  : %8.4f deg  @ t = %.4f s', v, t(j));
[v, j] = max(height);
R = pushLine(R, '  高度最大                : %8.4f m    @ t = %.4f s', v, t(j));
[v, j] = max(thrustWorst);
R = pushLine(R, '  推力占比最大            : %8.4f %%   @ t = %.4f s', v, t(j));
% ★ 张力是 n x nSteps 矩阵，必须用**线性下标**取全局最小值，再换算成 (绳号, 时刻)。
%   不能直接写 t(j)：j 是线性下标，最大可达 n*nSteps > nSteps ⇒ 越界报错。
[vMin, jLin] = min(T(:));
iLinkMin = mod(jLin - 1, n) + 1;
jStepMin = floor((jLin - 1) / n) + 1;
R = pushLine(R, '  张力最小                : %8.4f N    @ t = %.4f s（第 %d 根绳）', ...
    vMin, t(jStepMin), iLinkMin);
R = pushLine(R, '  各绳张力范围            :');
for i = 1:n
    R = pushLine(R, '     绳 %d: [%.6g, %.6g] N', i, min(T(i, :)), max(T(i, :)));
end
R = pushLine(R, '  稳态张力（末 20%% 绷紧段）: [%s] N', ...
    numVecText(mean(T(:, steadyTautIdx), 2).', '%.6g '));
R = pushLine(R, '  稳态期望张力            : [%s] N', ...
    numVecText(mean(desT(:, steadyTautIdx), 2).', '%.6g '));
R = pushLine(R, '  稳态绳向竖直流向分量    : [%s]（悬停应≈+1）', ...
    numVecText(mean(qVert(:, steadyTautIdx), 2).', '%.6g '));
R = pushLine(R, '  稳态绳向误差（ACTIVE）  : %.4f deg', mean(linkErrWorst(activeIdx)));
R = pushLine(R, '  稳态位置误差（ACTIVE）  : %.4f m', mean(errNorm(activeIdx)));
R = pushLine(R, '  稳态偏航误差（ACTIVE）  : %.4f deg', mean(abs(yawErrDeg(activeIdx))));
R = pushLine(R, '  稳态位置误差（末 20%%）  : %.4f m ← 含降落段，勿用', mean(errNorm(steadyIdx)));
R = pushLine(R, '  巡航段位置误差（ACTIVE 后 2/3）: %.4f m', ...
    mean(errNorm(activeIdx(max(1, round(numel(activeIdx) / 3)):end))));
R = pushLine(R, '  绳长范数最大偏差        : %.6g（||q_i|| 应恒为 1）', max(abs(qNorm(:) - 1)));
R = pushLine(R, '  状态码经历              : [%s]', numVecText(unique(modeSeq), '%d '));
if isfield(sim, 'nonFiniteSolveCount')
    R = pushLine(R, '  ★ 非有限解步数          : %d（必须为 0，唯一不会骗人的发散判据）', ...
        sim.nonFiniteSolveCount);
end
if isfield(sim, 'onsetStep') && ~isnan(sim.onsetStep)
    R = pushLine(R, '  ★ 发散起点步号          : %d（t = %.4f s）', sim.onsetStep, t(min(max(sim.onsetStep, 1), nSteps)));
end

R = pushLine(R, '');
R = pushLine(R, '【F】sim.summary 全字段');
R = pushLine(R, '  ⚠ minLinkVerticalComponent = 0 与"绳长范数偏差 = 1"来自 SLACK / TAKEUP /');
R = pushLine(R, '    LANDING_RELEASE 这些**绳索未建模**的相位（那时 q_i 记作 0 向量），');
R = pushLine(R, '    不是动力学问题 —— 绳向类判据必须限定在绷紧段（模式 2/3/4）。');
if isfield(sim, 'summary')
    R = dumpStructInto(R, sim.summary, 'summary');
else
    R = pushLine(R, '  （无 summary 字段）');
end

R = pushLine(R, '');
R = pushLine(R, '==============================================================================');

% ------------------------------------------------------------------ 落盘
reportPath = fullfile(outDir, 'dump_report.txt');
writeTextFile(reportPath, R);

writeTableCsv(fullfile(outDir, 'dump_trajectory.csv'), ...
    {'time', 'modeCode', 'tensionScale', ...
     'pos_x', 'pos_y', 'pos_z', 'height', ...
     'des_x', 'des_y', 'des_z', 'desHeight', ...
     'err_x', 'err_y', 'err_z', 'errNorm'}, ...
    [t; modeSeq; tScale; loadPos; height; desPos; desH; errVec; errNorm]);

writeTableCsv(fullfile(outDir, 'dump_velocity.csv'), ...
    {'time', 'vel_x', 'vel_y', 'vel_z', 'speed'}, ...
    [t; loadVel; sqrt(sum(loadVel .^ 2, 1))]);

tensionHeader = [{'time'}, linkNames('T', n), linkNames('desT', n), ...
                 linkNames('thrustPct', n)];
writeTableCsv(fullfile(outDir, 'dump_tension.csv'), tensionHeader, ...
    [t; T; desT; thrustPct]);

linkHeader = [{'time'}, linkNames('linkErrDeg', n), linkNames('ropeDist', n), ...
              linkNames('ropeSlack', n), linkNames('qVert', n)];
writeTableCsv(fullfile(outDir, 'dump_link.csv'), linkHeader, ...
    [t; linkErrDeg; sim.ropeDistanceLog; sim.ropeSlackLog; qVert]);

writeTableCsv(fullfile(outDir, 'dump_attitude.csv'), ...
    {'time', 'loadAttErr_x', 'loadAttErr_y', 'loadAttErr_z', ...
     'loadRate_x', 'loadRate_y', 'loadRate_z', ...
     'yawRef_deg', 'yaw_deg', 'yawErr_deg'}, ...
    [t; sim.loadAttitudeErrorLog; loadRate; yawRefDeg; yawDeg; yawErrDeg]);

vehHeader = [{'time'}, linkNames('vehX', n), linkNames('vehY', n), linkNames('vehZ', n)];
vehPos = sim.vehiclePositionLog;                     % 3 x n x nSteps
writeTableCsv(fullfile(outDir, 'dump_vehicle.csv'), vehHeader, ...
    [t; reshape(vehPos(1, :, :), n, nSteps); ...
        reshape(vehPos(2, :, :), n, nSteps); ...
        reshape(vehPos(3, :, :), n, nSteps)]);

% ------------------------------------------------------------------ 控制台
fprintf('\n');
for k = 1:numel(R)
    fprintf('%s\n', R{k});
end
fprintf('\n[dump] 完成，文件已写入：%s\n', outDir);
end

% ======================================================================
%                            局部函数
% ======================================================================

function R = pushLine(R, fmt, varargin)
% 追加一行到报告（cellstr）。
R{end + 1} = sprintf(fmt, varargin{:});
end

function R = dumpStructInto(R, s, prefix)
% 递归展开 struct（含标量 struct 数组 / 非标量 struct 数组），逐行 name = value。
if ~isstruct(s)
    R = pushLine(R, '%s = <%s>', prefix, class(s));
    return
end
if numel(s) > 1
    for k = 1:numel(s)
        R = dumpStructInto(R, s(k), sprintf('%s(%d)', prefix, k));
    end
    return
end
names = fieldnames(s);
for k = 1:numel(names)
    fullName = sprintf('%s.%s', prefix, names{k});
    v = s.(names{k});
    if isstruct(v)
        R = dumpStructInto(R, v, fullName);
    elseif ischar(v)
        R = pushLine(R, '  %-40s = %s', fullName, v);
    elseif isempty(v)
        R = pushLine(R, '  %-40s = []', fullName);
    elseif islogical(v) && isscalar(v)
        R = pushLine(R, '  %-40s = %d', fullName, v);
    elseif isnumeric(v) && isscalar(v)
        R = pushLine(R, '  %-40s = %.10g', fullName, v);
    elseif isnumeric(v)
        flat = v(:).';
        if numel(flat) > 12
            R = pushLine(R, '  %-40s = [%d x %d] 前 6 个: %s ...', ...
                fullName, size(v, 1), size(v, 2), numVecText(flat(1:6), '%.6g '));
        else
            R = pushLine(R, '  %-40s = [%s]', fullName, numVecText(flat, '%.6g '));
        end
    else
        R = pushLine(R, '  %-40s = <%s>', fullName, class(v));
    end
end
end

function [peak, peakTime] = peakInWindow(x, t, t0, t1)
% 行向量 x 在时间窗 [t0,t1] 内的最大值及其时刻。
idx = find(t >= t0 & t <= t1);
if isempty(idx)
    peak = NaN; peakTime = NaN;
    return
end
[peak, j] = max(x(idx));
peakTime = t(idx(j));
end

function [value, valueTime, linkIndex] = peakInWindowMatrix(X, t, t0, t1, how)
% X 为 n x nSteps：求窗口内的极值（how = 'min' 或 'max'）及其时刻与绳编号。
idx = find(t >= t0 & t <= t1);
if isempty(idx)
    value = NaN; valueTime = NaN; linkIndex = NaN;
    return
end
if strcmpi(how, 'min')
    [perStep, linkOfStep] = min(X(:, idx), [], 1);
    [value, j] = min(perStep);
else
    [perStep, linkOfStep] = max(X(:, idx), [], 1);
    [value, j] = max(perStep);
end
linkIndex = linkOfStep(j);
valueTime = t(idx(j));
end

function writeTableCsv(path, header, data)
% 写 CSV。data 必须是 M x N（M = 变量个数 = numel(header)，N = 采样数）。
% ★ fprintf 按**列优先**顺序取参数，所以 data 必须是 M x N：
%   这样遍历顺序正好是"第 1 个采样的 M 个变量、第 2 个采样…"，与格式串匹配。
if size(data, 1) ~= numel(header)
    error('dump_sim_data:BadTableShape', ...
        '%s：data 有 %d 行，而表头有 %d 列（应为 M x N，M = 变量数）。', ...
        path, size(data, 1), numel(header));
end
fid = fopen(path, 'w');
if fid < 0
    error('dump_sim_data:OpenFailed', '无法写入 %s', path);
end
fprintf(fid, '%s\n', strjoin(header, ','));
fmt = [repmat('%.10g,', 1, size(data, 1) - 1), '%.10g\n'];
fprintf(fid, fmt, data);
fclose(fid);
fprintf('[dump] 写入 %-28s (%d 变量 x %d 采样)\n', ...
    fileNameOf(path), size(data, 1), size(data, 2));
end

function writeTextFile(path, lines)
% 写文本报告。
fid = fopen(path, 'w');
if fid < 0
    error('dump_sim_data:OpenFailed', '无法写入 %s', path);
end
for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
end
fclose(fid);
fprintf('[dump] 写入 %-28s (%d 行)\n', fileNameOf(path), numel(lines));
end

function name = fileNameOf(path)
[~, base, ext] = fileparts(path);
name = [base, ext];
end

function names = linkNames(prefix, n)
% 生成 {'T1','T2',...} 这类表头单元数组。
names = cell(1, n);
for i = 1:n
    names{i} = sprintf('%s%d', prefix, i);
end
end

function out = numVecText(v, fmt)
% 把数值向量格式化成一行文本：逐项 sprintf(fmt) 后用空格连接。
% 相比内置的 num2str(向量, 格式串)：
%   ① 能保留 '%+.4f' 这类**显式正号**（num2str 不保证）；
%   ② 调用处习惯在格式串尾部多写一个空格，这里自动吃掉，避免"每项后缀一个空格"。
if ~isempty(fmt) && fmt(end) == ' '
    fmt(end) = [];
end
parts = cell(1, numel(v));
for i = 1:numel(v)
    parts{i} = sprintf(fmt, v(i));
end
out = strjoin(parts, ' ');
end

function v = getOr(s, name, fallback)
if isfield(s, name)
    v = s.(name);
else
    v = fallback;
end
end
