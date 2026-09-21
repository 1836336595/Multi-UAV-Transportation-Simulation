function crazyflie_slung_diagnose()
%CRAZYFLIE_SLUNG_DIAGNOSE  发散定位脚本（不画图，只报告）
%
% 用途：当 `crazyflie_slung_demo` 出现
%       「矩阵为奇异值…RCOND = NaN」这类警告（= 状态里已进入 NaN）
%       或图形阶段报"范围必须为包含递增的数值的 2 元素向量"时，
%       用本脚本一次性给出：
%         1) **实际生效的关键参数** —— 用来排除"MATLAB 路径上是旧版参数文件"
%            这一类问题（参数变了但报错症状没变，往往就是这个原因）；
%         2) **第一个非有限步的位置与时刻**；
%         3) 该步**前后各通道的数值**（负载位置/角速度、绳向、机体角速度、
%            推力、张力），用于判断是哪一路先炸的。
%
% 用法（在 MATLAB 命令行）：
%     >> cd('...\v3')          % 进到本工程目录
%     >> crazyflie_slung_diagnose
%
% 说明：本脚本只读不写，不会覆盖任何文件。

fprintf('==============================================================\n');
fprintf(' 发散定位诊断\n');
fprintf('==============================================================\n');

% ------------------------------------------------ 1) 实际生效的参数
cfg = crazyflie_slung_parameters();
fe = cfg.figureEight;

fprintf('\n[1] 实际生效的关键参数（请对照 README §11 核对）\n');
fprintf('  ---- 挂点 attachPoints（列 = 挂点 rho_i）----\n');
disp(cfg.payload.attachPoints);
fprintf('  cond(P*P'') 相关：挂点质心 = (%.4f, %.4f)\n', ...
    mean(cfg.payload.attachPoints(1, :)), mean(cfg.payload.attachPoints(2, :)));

fprintf('  ---- 机体姿态环 ----\n');
fprintf('  attitudeController.kR    = [%g; %g; %g]\n', cfg.attitudeController.kR);
fprintf('  attitudeController.kOmega= [%g; %g; %g]\n', cfg.attitudeController.kOmega);
fprintf('  rateLoop.bandwidth       = [%g; %g; %g]\n', cfg.rateLoop.bandwidth);
fprintf('  等效姿态环 wn = sqrt(BW*kR) = %.1f rad/s, zeta = %.3f\n', ...
    sqrt(cfg.rateLoop.bandwidth(1) * cfg.attitudeController.kR(1)), ...
    cfg.rateLoop.bandwidth(1) * (1 + cfg.attitudeController.kOmega(1)) ...
    / (2 * sqrt(cfg.rateLoop.bandwidth(1) * cfg.attitudeController.kR(1))));

fprintf('  ---- 八字工况 ----\n');
fprintf('  amplitudeX=%.2f  amplitudeY=%.2f  cycles=%.1f\n', ...
    fe.amplitudeX, fe.amplitudeY, fe.cycles);
fprintf('  takeoff=%.1f  cruise=%.1f  landing=%.1f  duration=%.1f\n', ...
    fe.takeoffDuration, fe.cruiseDuration, fe.landingDuration, cfg.simulation.duration);
fprintf('  blendTime=%.2f  cruiseHeight=%.2f  lockYaw=%d\n', ...
    fe.blendTime, fe.cruiseHeight, fe.lockYaw);

fprintf('  ---- 负载控控制/绳向控器 ----\n');
fprintf('  kx=[%s]  kv=[%s]  ki=[%s]\n', mat2str(cfg.loadController.kx.', 4), ...
    mat2str(cfg.loadController.kv.', 4), mat2str(cfg.loadController.ki.', 4));
fprintf('  kq=%.1f  komega=%.1f  kqIntegral=%.1f\n', ...
    cfg.linkController.kq, cfg.linkController.komega, cfg.linkController.kqIntegral);

fprintf('  ---- 初始条件 ----\n');
fprintf('  initial.position   = [%s]\n', mat2str(cfg.initial.position.', 4));
fprintf('  initial.thrustNewton = [%s]   (1 x n = %d)\n', ...
    mat2str(cfg.initial.thrustNewton(:).', 6), numel(cfg.initial.thrustNewton));
fprintf('  initial.linkUnits  =\n'); disp(cfg.initial.linkUnits);

fprintf('  ---- 障碍物 ----\n');
fprintf('  obstacles.positions =\n'); disp(cfg.obstacles.positions);
fprintf('  height=%.2f  radius=%.2f  baseZ=%.2f\n', ...
    cfg.obstacles.height, cfg.obstacles.radius, cfg.obstacles.baseZ);

% ------------------------------------------------ 2) 跑仿真（不画图）
fprintf('\n[2] 运行仿真（不绘图）…\n');

% ★★ 必须清掉 crazyflie_slung_dynamics 里的 persistent 标志！
%   `crazyflie_slung_dynamics` 里有个 persistent 变量，用来保证"非有限量详细报告"
%   只打印一次（否则 26000 步会刷屏）。但 persistent 在**同一个 MATLAB 会话内不会
%   自动重置** —— 于是同一个诊断脚本连跑两次，第二次会**静默**，
%   让人误判成"这次没发散"。这个坑我本人踩过（两次运行结果不一致）。
% ★ 用函数形式 clear('...')：`clear crazyflie_slung_dynamics` 那种命令式写法
%   会把函数名当成变量名，静态检查的"先用后定义"会误报。
clear('crazyflie_slung_dynamics')

cfg.visualization.plot = false;
cfg.visualization.animate = false;
cfg.visualization.saveVideo = false;

t0 = tic;
try
    sim = crazyflie_slung_simulation(cfg);
catch err
    fprintf('  仿真本身抛出异常：\n');
    fprintf('    %s\n', err.message);
    if ~isempty(err.stack)
        fprintf('    位置：%s (第 %d 行)\n', err.stack(1).name, err.stack(1).line);
    end
    return;
end
fprintf('  完成，用时 %.1f s\n', toc(t0));

% ------------------------------------------------ 3) 找第一个非有限步
logNames = {'loadPositionLog', 'loadVelocityLog', 'loadRotationLog', 'loadBodyRateLog', ...
    'linkUnitLog', 'linkRateLog', 'rotationLog', 'bodyRateLog', ...
    'thrustLog', 'tensionLog', 'momentLog', 'positionErrorLog'};

nSteps = numel(sim.time);
bad = false(1, nSteps);
firstBadByLog = struct();
for j = 1:numel(logNames)
    nm = logNames{j};
    if ~isfield(sim, nm) || isempty(sim.(nm))
        continue;
    end
    v = sim.(nm);
    nLast = size(v, ndims(v));
    if nLast ~= nSteps
        continue;                                  % 形状不符，跳过
    end
    m = reshape(v, [], nLast);
    badHere = any(~isfinite(m), 1);
    bad = bad | badHere;
    k = find(badHere, 1, 'first');
    if ~isempty(k)
        firstBadByLog.(nm) = k;
    end
end

fprintf('\n[3] 非有限状态定位\n');

% ★★ 先看"坏步计数"。这是**唯一不会骗人的判据**：
%   如果 6x6 代数系统曾出现非有限量，坏步的解会被置零，NaN 不会传进状态，
%   于是**所有状态日志仍然全是有限值** —— 只看日志会把"发散"误判成"正常"。
nBadStep = nan;      % 用小写 nan（与 whitelist 一致；NaN 同义）
if isfield(sim, 'summary') && isfield(sim.summary, 'nonFiniteSolveCount')
    nBadStep = sim.summary.nonFiniteSolveCount;
end
if isfinite(nBadStep)
    if nBadStep == 0
        fprintf('  ✅ 坏步计数 = 0：6x6 代数系统全程有限，闭环真的没发散。\n');
    else
        fprintf('  ❌ 坏步计数 = %d 步！代数系统曾出现非有限量（见上方一次性报告）。\n', nBadStep);
        fprintf('     注意：这类坏步的解被置零保护，**状态日志可能仍然全部有限**，\n');
        fprintf('     所以不能只看日志。病因已由 crazyflie_slung_dynamics 打印。\n');
    end
else
    fprintf('  （本版 sim.summary 没有 nonFiniteSolveCount 字段，跳过；旧版无此保护）\n');
end
% ★ 发散**起点**步号（循环内的"起点探测"会打印该拍的上文）
if isfield(sim, 'summary') && isfield(sim.summary, 'onsetStep') ...
        && isfinite(sim.summary.onsetStep)
    k0 = sim.summary.onsetStep;
    fprintf('  ★ 发散起点：k = %d, t = %.4f s（该拍的上文已在 [2] 段打印）\n', ...
        k0, sim.time(min(max(k0, 1), numel(sim.time))));
    fprintf('     换算：k·dt = %.4f s；T1 = 3.0 s、环绕段 [3, 48)。\n', k0 * cfg.simulation.dt);
end

% ★★ 时间历程抽样 —— 用来判"从哪一拍开始跑掉"。
%   比设阈值更可靠：阈值总会漏（实测把阈值设成 ||Omega0||>1e5 时，
%   负载已经飞到 145 m 外了，而更早的"位置跑掉"完全没被触发）。
printTimeHistory(sim, cfg);

if ~any(bad)
    fprintf('  日志有限性：全程有限。\n');
    if isfinite(nBadStep) && nBadStep == 0
        fprintf('  ⇒ 本次运行**确实没有发散**。若 demo 仍报 RCOND = NaN，\n');
        fprintf('    请确认运行的是本目录下的 .m 文件（which -all crazyflie_slung_parameters）。\n');
    end
    return;
end

kBad = find(bad, 1, 'first');
fprintf('  ❌ 第一个非有限步：k = %d，t = %.3f s\n', kBad, sim.time(kBad));
fprintf('     各日志首次出问题的步：\n');
fn = fieldnames(firstBadByLog);
for j = 1:numel(fn)
    fprintf('       %-22s  k = %d  (t = %.3f s)\n', fn{j}, ...
        firstBadByLog.(fn{j}), sim.time(firstBadByLog.(fn{j})));
end

% 该步前后各通道数值
kShow = max(1, kBad - 1);
fprintf('\n[4] 故障前一步 (k = %d, t = %.3f s) 的数值\n', kShow, sim.time(kShow));
showState(sim, cfg, kShow);
fprintf('\n[5] 故障当步 (k = %d, t = %.3f s) 的数值\n', kBad, sim.time(kBad));
showState(sim, cfg, kBad);

fprintf('\n提示：把上面 [1] ~ [5] 的完整输出回贴，即可定位是哪一路先炸的。\n');
end

% ======================================================================
function showState(sim, cfg, k)
n = cfg.vehicle.count;
fprintf('   负载位置      : [%s]\n', mat2str(sim.loadPositionLog(:, k).', 5));
fprintf('   负载速度      : [%s]\n', mat2str(sim.loadVelocityLog(:, k).', 5));
fprintf('   位置误差      : %.6f m\n', sim.positionErrorLog(k));
fprintf('   负载姿态 R0   :\n'); disp(sim.loadRotationLog(:, :, k));
fprintf('   负载体角速度  : [%s] rad/s\n', mat2str(sim.loadBodyRateLog(:, k).', 5));
fprintf('   绳向 q_i      :\n'); disp(sim.linkUnitLog(:, :, k));
fprintf('   绳向角速率    : [%s] rad/s\n', mat2str(sim.linkRateLog(:, k).', 5));
fprintf('   ||e_q|| (deg) : [%s]\n', mat2str(sim.linkErrorLog(:, k).', 5));
fprintf('   机体角速度    : [%s] rad/s\n', mat2str(sim.bodyRateLog(:, k).', 5));
fprintf('   角速度指令    : [%s] rad/s\n', mat2str(sim.omegaCommandLog(:, k).', 5));
fprintf('   实际推力      : [%s] N  (上限 %.4f)\n', ...
    mat2str(sim.thrustLog(:, k).', 5), cfg.vehicle.maxTotalThrust);
fprintf('   推力占比      : [%s] %%\n', mat2str(sim.thrustPctLog(:, k).', 4));
fprintf('   绳索张力      : [%s] N\n', mat2str(sim.tensionLog(:, k).', 5));
fprintf('   （n = %d）\n', n);
end

% ======================================================================
function printTimeHistory(sim, cfg)
% 抽样打印时间历程，用于定位"从哪一拍开始跑掉"。
% ★ 为什么要有这个：单靠"阈值探测"一定会漏。实测把阈值设成 ||Omega0|| > 1e5，
%   触发时负载已经飞到 145 m 外 —— 更早的"位置跑掉"完全没被抓住。
%   时间历程是**无阈值**的证据：直接看哪一拍开始偏离。
nSteps = numel(sim.time);
fprintf('\n[3.5] 时间历程抽样（定位"从哪一拍开始跑掉"）\n');
fprintf('           k       t(s)      |e| (m)    |Om0|(rad/s)   Om0_z(rad/s)   maxThr%%     max|T|(N)\n');
fprintf('    ---- 全程粗扫 ----\n');
printHistoryRows(sim, unique(round(linspace(1, nSteps, 36))));
nFine = max(2, round(0.3 * nSteps));
fprintf('    ---- 前 30%%（0 ~ %.2f s）细扫 ----\n', sim.time(nFine));
printHistoryRows(sim, unique(round(linspace(1, nFine, 36))));
end

% ======================================================================
function printHistoryRows(sim, idx)
% ★ 一定要打印 Om0_z：本次故障的**真正病根**就是偏航转速失控
%   （实测起点那拍 Omega0 = [-0.0009, 0.023, -5.01]，只有 z 分量炸了），
%   而只看 ||Omega0|| 会被 x/y 的小量稀释、看不出方向。
for kk = idx
    e = nan; o = nan; th = nan; te = nan; oz = nan;
    if isfield(sim, 'positionErrorLog') && numel(sim.positionErrorLog) >= kk
        e = sim.positionErrorLog(kk);
    end
    if isfield(sim, 'loadBodyRateLog') && size(sim.loadBodyRateLog, 2) >= kk
        o = norm(sim.loadBodyRateLog(:, kk));
        oz = sim.loadBodyRateLog(3, kk);
    end
    if isfield(sim, 'thrustPctLog') && size(sim.thrustPctLog, 2) >= kk
        th = max(sim.thrustPctLog(:, kk));
    end
    if isfield(sim, 'tensionLog') && size(sim.tensionLog, 2) >= kk
        te = max(sim.tensionLog(:, kk));
    end
    fprintf('    %7d  %9.4f  %11.4e  %12.4e  %11.4e  %9.1f  %12.4e\n', ...
        kk, sim.time(kk), e, o, oz, th, te);
end
end
