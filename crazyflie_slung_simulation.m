function sim = crazyflie_slung_simulation(userCfg)
%CRAZYFLIE_SLUNG_SIMULATION n 架 Crazyflie 2.1 Brushless 协同吊运刚体负载仿真。
%
% 用法：
%   sim = crazyflie_slung_simulation();                 % 默认：2 组曲线 + 三维动画
%   cfg.visualization.animate = false;
%   sim = crazyflie_slung_simulation(cfg);              % 只看曲线
%   sim = crazyflie_slung_simulation(struct());         % 全默认
%
% 仿真结构：
%   [负载 x0, R0] + n 组 [绳索 q_i] + [机体 R_i]
%            |
%   负载外环 (20)-(21) -> Fd, Md
%            |
%   张力分配 (13)(22)-(23) -> mu_id（n 根绳索的期望张力）
%            |
%   绳向 (24)(25) + 平行分量 (17) + 绳向环 (27) -> u_i
%            |
%   推力 f_i = -u_i' R_i e3  +  姿态外环 -> Omega_cmd_i
%            |
%   推力一阶执行器 + 角速度内环（等效固件速率环）-> 力矩 M_i
%            |
%   完整动力学 (5)-(8) 积分
%
% 本文件只做：初始化、调用控制器、模拟内外环执行机构、积分动力学、记录数据。

if nargin < 1 || isempty(userCfg)
    userCfg = struct();
end

cfg = crazyflie_slung_parameters(userCfg);
dt = cfg.simulation.dt;
nSteps = floor(cfg.simulation.duration / dt) + 1;
time = (0:nSteps - 1) * dt;
n = cfg.vehicle.count;
rhoAll = cfg.payload.attachPoints;      % 3 x n

% ------------------------------------------------------------------ 初始状态
state = struct();
state.loadPosition = cfg.initial.position(:);
state.loadVelocity = cfg.initial.velocity(:);
state.loadRotation = projectSO3(cfg.initial.R0);
state.loadBodyRate = cfg.initial.loadBodyRate(:);
state.linkUnits = zeros(3, n);
state.linkRates = zeros(3, n);
for i = 1:n
    state.linkUnits(:, i) = normalizeVector(cfg.initial.linkUnits(:, i));
    state.linkRates(:, i) = cfg.initial.linkRates(:, i);
end
state.rotations = zeros(3, 3, n);
state.bodyRates = zeros(3, n);
state.bodyTorques = zeros(3, n);
for i = 1:n
    state.rotations(:, :, i) = projectSO3(cfg.initial.vehicleR(:, :, i));
    state.bodyRates(:, i) = cfg.initial.bodyRates(:, i);
end
state.loadAcceleration = zeros(3, 1);
state.loadBodyAcceleration = zeros(3, 1);

% 各机当前实际推力（一阶执行器状态）
% ★ thrustNewton 允许是标量（旧写法，各机相同）或 1 x n 向量（按挂点几何解出的
%   不等悬停推力）。挂点不对称时三根绳张力是 2:1:1，必须用向量给定，
%   否则起步瞬间会有一次 100% 推力饱和（实测）。
if isscalar(cfg.initial.thrustNewton)
    actualThrust = cfg.initial.thrustNewton * ones(1, n);
else
    actualThrust = reshape(cfg.initial.thrustNewton, 1, n);
end

% 控制器跨步状态
memory = struct();
memory.positionIntegral = zeros(3, 1);
memory.linkIntegrals = zeros(3, n);
memory.previousLinkUnits = [];

% ------------------------------------------------------------------ 日志分配
sim = struct();
sim.time = time;
sim.config = cfg;
sim.loadPositionLog = zeros(3, nSteps);
sim.loadVelocityLog = zeros(3, nSteps);   % ★ 状态量，必须记录（见下面的说明）
sim.loadRotationLog = zeros(3, 3, nSteps);
sim.linkUnitLog = zeros(3, n, nSteps);
sim.vehiclePositionLog = zeros(3, n, nSteps);
sim.rotationLog = zeros(3, 3, n, nSteps);
sim.bodyRateLog = zeros(3, n, nSteps);
sim.linkRateLog = zeros(3, n, nSteps);
sim.parallelForceLog = zeros(3, n, nSteps);
sim.perpendicularForceLog = zeros(3, n, nSteps);
sim.totalForceLog = zeros(3, n, nSteps);
sim.commandLog = repmat(emptyCommandLog(n), 1, nSteps);
sim.tensionLog = zeros(n, nSteps);
sim.thrustLog = zeros(n, nSteps);
sim.thrustPctLog = zeros(n, nSteps);
sim.momentLog = zeros(3, n, nSteps);
sim.attitudeErrorLog = zeros(n, nSteps);
sim.linkErrorLog = zeros(n, nSteps);
sim.positionErrorLog = zeros(1, nSteps);
sim.positionErrorVectorLog = zeros(3, nSteps);
sim.omegaCommandLog = zeros(3, n, nSteps);
sim.desiredTensionLog = zeros(3, n, nSteps);
% 负载姿态误差的按轴日志（3 x nSteps）：用于区分可控轴（roll/pitch）
% 与欠驱动轴（yaw，见 README §5.1）。负载只有一个，故不重复 n 份。
sim.loadAttitudeErrorLog = zeros(3, nSteps);
% 负载角速度日志（3 x nSteps）：用于检查偏航漂移率。
sim.loadBodyRateLog = zeros(3, nSteps);
% ★★★ 期望负载 yaw 日志（1 x nSteps，单位 rad）：取**参考姿态 R0d 的第一轴方位角**。
%   必须在这里记录，不能在可视化里事后重调 cfg.referenceFcn —— 
%   匿名句柄的**多输出**调用在部分 MATLAB 版本上会失败（本文件顶部的
%   callReferenceFunction 之所以要 6→5→4→3 逐级回退就是为了这个），
%   事后调用一旦失败就会静默退化成"只画实际 yaw"，用户看不到参考曲线。
%   在仿真循环里记录则走的是同一条已验证的分派路径，不会失败。
sim.loadYawRefLog = zeros(1, nSteps);

previousLoadAcceleration = zeros(3, 1);
previousLoadBodyAcceleration = zeros(3, 1);
% ★ "6x6 代数系统出现非有限量"的步数累计（见动力学里 info.nonFiniteInputs 的说明）。
%   这是唯一不会骗人的发散判据：坏步的解被置零后 NaN 不进状态，
%   状态日志可能全有限，只看日志会误判"正常"。
nonFiniteStepCount = 0;
% ★ 发散起点探测（见循环内说明）：一次性的起点报告标志与起点步号
onsetReported = false;
onsetStep = nan;    % 用小写 nan（与 whitelist 一致）

for k = 1:nSteps
    t = time(k);
    desired = referenceState(t, cfg);

    % -------- 控制器 --------
    [command, memory] = crazyflie_slung_controller(state, desired, memory, cfg);
    sim.commandLog(k) = logCommand(command, n);
    sim.tensionLog(:, k) = command.tensions(:);
    sim.linkErrorLog(:, k) = sqrt(sum(command.linkDirectionErrors.^2, 1)).';
    sim.attitudeErrorLog(:, k) = sqrt(sum(command.attitudeErrors.^2, 1)).';
    % 负载姿态误差按轴记录。command.loadAttitudeError 已由控制器给出
    % （论文 (21) 的 e_R0 = 0.5 vee(R0d' R0 - R0' R0d)），是 3x1 向量。
    if isfield(command, 'loadAttitudeError') && numel(command.loadAttitudeError) == 3
        sim.loadAttitudeErrorLog(:, k) = command.loadAttitudeError(:);
    end
    sim.loadBodyRateLog(:, k) = state.loadBodyRate;
    % 期望负载 yaw = R0d 第一轴的方位角（本工况下 = 参考运动方向）
    sim.loadYawRefLog(k) = atan2(desired.rotation(2, 1), desired.rotation(1, 1));
    sim.positionErrorLog(k) = norm(command.positionError);
    sim.positionErrorVectorLog(:, k) = command.positionError;
    sim.omegaCommandLog(:, :, k) = command.omegaCommands;
    sim.desiredTensionLog(:, :, k) = command.desiredTensions;

    % 记录当前状态
    sim.loadPositionLog(:, k) = state.loadPosition;
    % ★ loadVelocity 是**状态量**，必须记录：它参与积分（loadPosition 由它积分而来），
    %   漏记会让"日志里有没有 NaN"这类全局自检**看不到它** ——
    %   排障时曾因此把"状态已坏但日志全有限"误判为"没有发散"。
    sim.loadVelocityLog(:, k) = state.loadVelocity;
    sim.loadRotationLog(:, :, k) = state.loadRotation;
    sim.linkUnitLog(:, :, k) = state.linkUnits;
    sim.linkRateLog(:, :, k) = state.linkRates;
    sim.rotationLog(:, :, :, k) = state.rotations;
    sim.bodyRateLog(:, :, k) = state.bodyRates;
    sim.parallelForceLog(:, :, k) = command.parallelForces;
    sim.perpendicularForceLog(:, :, k) = command.perpendicularForces;
    sim.totalForceLog(:, :, k) = command.totalForces;
    sim.thrustPctLog(:, k) = command.thrustPercentage(:);
    for i = 1:n
        sim.vehiclePositionLog(:, i, k) = state.loadPosition ...
            + state.loadRotation * rhoAll(:, i) ...
            - cfg.link.length * state.linkUnits(:, i);
    end

    if k == nSteps
        break;
    end

    % -------- 推力执行器：一阶响应（逐机独立） --------
    alpha = min(1, dt / max(cfg.vehicle.thrustTimeConstant, eps));
    actualThrust = actualThrust + alpha * (command.thrustDesired(:).' - actualThrust);
    actualThrust = min(max(actualThrust, 0), cfg.vehicle.maxTotalThrust);
    sim.thrustLog(:, k + 1) = actualThrust(:);

    % -------- 角速度内环：等效 Crazyflie 固件速率环（逐机独立） --------
    % 一阶闭环等效：Omega_dot_i = K_i*(Omega_cmd_i - Omega_i)，
    % 对应力矩 M_i = J_i*Omega_dot_i + Omega_i x J_i Omega_i
    for i = 1:n
        rateError = command.omegaCommands(:, i) - state.bodyRates(:, i);
        bodyRateDotCmd = cfg.rateLoop.bandwidth .* rateError;
        Omegai = state.bodyRates(:, i);
        moment = cfg.vehicle.inertia * bodyRateDotCmd ...
            + cross(Omegai, cfg.vehicle.inertia * Omegai);
        state.bodyTorques(:, i) = moment;
        sim.momentLog(:, i, k + 1) = moment;
    end

    % -------- 实际作用力：物理模型 -f_i R_i e3 --------
    % 真实四旋翼的合力方向由**实际机体姿态 R_i** 决定，幅值由实际推力决定。
    % 只有当姿态收敛到期望姿态 R_ic 时才有 -f_i R_i e3 -> u_i。
    % 不要写成 u_i * (f_actual/f_cmd)：那会把力锁在理想方向，架空姿态环，
    % 且当机体大倾斜时 f_cmd -> 0 会把力幅值放大到无穷。
    uActual = zeros(3, n);
    for i = 1:n
        uActual(:, i) = -actualThrust(i) * (state.rotations(:, :, i) * [0; 0; 1]);
    end

    % -------- 动力学（论文 (5)-(8)） --------
    state.loadAcceleration = previousLoadAcceleration;
    state.loadBodyAcceleration = previousLoadBodyAcceleration;
    [derivative, dynInfo] = crazyflie_slung_dynamics(state, uActual, cfg);
    % ★ 累计"6x6 代数系统出现非有限量"的步数。这是**唯一不会骗人的发散判据**：
    %   坏步的解会被置零，NaN 不会传进状态，所以状态日志可能全部有限 ——
    %   只看 sim.*Log 里有没有 NaN 会把"发散"误判成"正常"（已实测踩过）。
    nonFiniteStepCount = nonFiniteStepCount + dynInfo.nonFiniteInputs;

    % ★★★ 发散"起点"探测（一次性报告）
    %   已确认的因果链：控制器输出 NaN ⇒ thrustDesired = NaN
    %   ⇒ `min(max(NaN,0), maxTotalThrust)` **把 NaN 洗成 0**（MATLAB 的 max/min 忽略 NaN）
    %   ⇒ actualThrust = 0 ⇒ uActual = 0 ⇒ 四旋翼完全不出力 ⇒ 负载失稳、翻转
    %   ⇒ 负载角速度爆到 1e289 ⇒ hat(Omega0)^2 平方溢出为 Inf ⇒ rhsBlock = NaN。
    %   所以必须先抓住"ThrustDesired 第一次非有限"的那一拍（= 真正的起点），
    %   而后面的坏步计数只是它的后果。
    if ~onsetReported
        onsetKind = '';
        onsetDetail = '';
        if ~all(isfinite(command.thrustDesired(:)))
            onsetKind = '控制指令 command.thrustDesired 出现非有限（NaN 起点在控制器）';
        elseif ~all(isfinite(command.desiredForce(:)))
            onsetKind = '期望合力 command.desiredForce 出现非有限';
        elseif ~all(isfinite(command.totalForces(:)))
            onsetKind = '控制力 command.totalForces (u_i) 出现非有限';
        elseif max(abs(actualThrust)) == 0
            onsetKind = '实际推力全为 0（指令已被 min/max 洗成 0）';
        elseif norm(state.loadBodyRate) > 5
            % ★ 阈值必须**很紧**：镜像全程 max|Omega0| 只有 0.86 rad/s，
            %   而可用的角速度指令上限是 5 rad/s。用 1e5 会漏掉真正的中早期发散
            %   （实测：负载已经飞到 145 m 外时 ||Omega0|| 仍 < 1e5）。
            onsetKind = '负载角速度 ||Omega0|| > 5 rad/s（远超正常 0.86）';
        elseif abs(state.loadPosition(3)) > 3 || norm(state.loadPosition(1:2)) > 5
            onsetKind = '负载跑出工作区（|z| > 3 m 或水平 > 5 m）';
        end
        if ~isempty(onsetKind)
            onsetReported = true;
            onsetStep = k;
            fprintf('\n*** 发散起点定位：k = %d, t = %.4f s ***\n', k, time(k));
            fprintf('    现象：%s\n', onsetKind);
            fprintf('    负载角速度 Omega0        = [%s]\n', mat2str(state.loadBodyRate.', 6));
            fprintf('    上一拍 Omega0_dot        = [%s]\n', mat2str(previousLoadBodyAcceleration.', 6));
            fprintf('    负载位置                 = [%s]\n', mat2str(state.loadPosition.', 6));
            fprintf('    期望推力 thrustDesired   = [%s]\n', mat2str(command.thrustDesired(:).', 6));
            fprintf('    实际推力 actualThrust    = [%s]\n', mat2str(actualThrust(:).', 6));
            fprintf('    期望合力 Fd              = [%s]\n', mat2str(command.desiredForce(:).', 6));
            fprintf('    期望力矩 Md              = [%s]\n', mat2str(command.desiredMoment(:).', 6));
            fprintf('    控制力 u_i (3x%d)        = [%s]\n', n, ...
                mat2str(reshape(command.totalForces, 1, []), 5));
            fprintf('    期望张力 mu_id (3x%d)    = [%s]\n', n, ...
                mat2str(reshape(command.desiredTensions, 1, []), 5));
            fprintf('    绳索张力                 = [%s]\n', mat2str(command.tensions(:).', 6));
            fprintf('    角速度指令               = [%s]\n', ...
                mat2str(reshape(command.omegaCommands, 1, []), 5));
            fprintf('    绳向 q_i (3x%d)          = [%s]\n', n, ...
                mat2str(reshape(state.linkUnits, 1, []), 5));
            fprintf('    ------------------------------------------------\n');
        end
    end

    % -------- 半隐式欧拉积分 --------
    state.loadVelocity = state.loadVelocity + dt * derivative.loadVelocity;
    state.loadPosition = state.loadPosition + dt * state.loadVelocity;
    state.loadBodyRate = state.loadBodyRate + dt * derivative.loadBodyRate;
    state.loadRotation = projectSO3(state.loadRotation ...
        * expSO3(state.loadBodyRate * dt));

    for i = 1:n
        state.linkRates(:, i) = state.linkRates(:, i) + dt * derivative.linkRates(:, i);
        state.linkUnits(:, i) = normalizeVector( ...
            state.linkUnits(:, i) + dt * state.linkRates(:, i));
        state.bodyRates(:, i) = state.bodyRates(:, i) + dt * derivative.bodyRates(:, i);
        state.bodyRates(:, i) = clampVector(state.bodyRates(:, i), ...
            -cfg.simulation.maxBodyRate, cfg.simulation.maxBodyRate);
        state.rotations(:, :, i) = projectSO3(state.rotations(:, :, i) ...
            * expSO3(state.bodyRates(:, i) * dt));

        % 数值安全：绳索角速度限幅
        rateNorm = norm(state.linkRates(:, i));
        if rateNorm > cfg.simulation.maxPendulumRate
            state.linkRates(:, i) = state.linkRates(:, i) ...
                * (cfg.simulation.maxPendulumRate / rateNorm);
        end
    end

    % -------- 绳索约束投影：保持 q_i' q_dot_i ≡ 0 --------
    % ★ 与"绳 vs 刚性连杆"直接相关，不要删 ★
    %
    % 绳在**绷紧**时与刚性连杆力学完全相同：都只沿 q_i 传轴向力。
    % 所以论文的动力学 (5)-(8)、张力分配 (13)(22)(23)、控制器 (27)(36)-(40)
    % 一个都不用改。唯一差别是绳多一条**单边约束** mu_i >= 0（只能拉不能推）。
    % 实测最小张力 0.1809 N（= 悬停张力 0.2616 N 的 69.2 %），全程未被激活，
    % 因此当前结果对一个绷紧的绳是物理正确的（见 README §9.7）。
    %
    % 本仿真里绳长不变量 |x_veh,i - 挂点| == l 是**构造性成立**的：
    % 无人机位置本身就用 veh = x0 + R0 rho_i - l q_i 定义（见下方记录段），
    % 实测偏差 2.2e-16 m（浮点精度）。所以绳长不需要额外投影。
    %
    % 这里唯一要做的是把 q_dot 的**径向分量**剥掉，使 q_i' q_dot_i ≡ 0。
    % 这是 ‖q_i‖ ≡ 1 求导的直接推论，也是论文 (1) q_dot = omega x q 的结论。
    % 旧代码只对 q 做 renormalize 而不投影 q_dot，会有 O(dt·|q_dot|^2) 的
    % 径向残留在每步注入噪声（实测 q'·q_dot 残差达 -2.1e-4）。
    for i = 1:n
        qi = state.linkUnits(:, i);
        state.linkRates(:, i) = state.linkRates(:, i) ...
            - qi * dot(qi, state.linkRates(:, i));
    end

    previousLoadAcceleration = derivative.loadAcceleration;
    previousLoadBodyAcceleration = derivative.loadBodyAcceleration;
end

% 补最后一步的记录
sim.thrustLog(:, nSteps) = actualThrust(:);
for i = 1:n
    sim.vehiclePositionLog(:, i, nSteps) = state.loadPosition ...
        + state.loadRotation * rhoAll(:, i) ...
        - cfg.link.length * state.linkUnits(:, i);
end

% ★ 把"坏步计数"挂到 sim 上再进 computeSummary。
%   ⚠ 不能直接写 `summary.x = nonFiniteStepCount` —— computeSummary 是**另一个函数**，
%     它看不到主函数工作区里的 nonFiniteStepCount（MATLAB 函数作用域是分开的）。
%     这正是静态检查第 5 类（先用后定义）能抓到的那种错。
sim.nonFiniteSolveCount = nonFiniteStepCount;
sim.onsetStep = onsetStep;
sim.summary = computeSummary(sim, cfg);

% 可视化（与数值计算解耦）
if cfg.visualization.plot || cfg.visualization.animate
    crazyflie_slung_visualization(sim, cfg);
end
end

% ======================================================================
function desired = referenceState(t, cfg)
% 读取静态目标或用户提供的轨迹函数。
if isempty(cfg.referenceFcn)
    desired.position = cfg.target.position;
    desired.velocity = cfg.target.velocity;
    desired.acceleration = cfg.target.acceleration;
    desired.rotation = cfg.target.R0;
    desired.bodyRate = cfg.reference.omegaD;
    desired.bodyRateDot = cfg.reference.omegaDotD;
else
    [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
        callReferenceFunction(cfg.referenceFcn, t);
    desired.position = position;
    desired.velocity = velocity;
    desired.acceleration = acceleration;
    desired.rotation = rotation;
    desired.bodyRate = bodyRate;
    desired.bodyRateDot = bodyRateDot;
end
desired.position = desired.position(:);
desired.velocity = desired.velocity(:);
desired.acceleration = desired.acceleration(:);
desired.rotation = projectSO3(desired.rotation);
desired.bodyRate = desired.bodyRate(:);
desired.bodyRateDot = desired.bodyRateDot(:);
end

function [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    callReferenceFunction(referenceFcn, t)
% 兼容三种轨迹函数契约：
%   ① 数值多输出（6 / 5 / 4 / 3 个输出）
%   ② 单输出 struct（含 position / velocity / acceleration / rotation 等字段）
%   ③ 单输出位置向量（退化契约）
%
% ★★ 陷阱（曾经报"此类型的变量不支持使用点进行索引"）★★
%   MATLAB 对**匿名函数句柄**执行 nargout() 恒返回 -1，因为输出个数"未知"。
%   而 cfg.referenceFcn 正是匿名句柄 @(t) crazyflie_slung_reference(t, cfg)，
%   所以实际走的是 nargout < 0 这条路径。旧版分派只写了 >= 6 / == 5 / == 4，
%   -1 三个分支都不命中 → 掉进最后的 struct 分支 → 把数值向量当结构体取字段
%   → 直接报错。（静态悬停工况走 cfg.referenceFcn = [] 的分支，
%   所以这条路径此前从未被执行，一直潜伏到现在。）
%
%   ⇒ 必须**显式**处理 nargout < 0，见下面的 ① 分支。

numberOfOutputs = nargout(referenceFcn);

if numberOfOutputs < 0
    % ① 匿名句柄：nargout 恒为 -1，只能逐级回退试探
    [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
        callReferenceAnonymous(referenceFcn, t);
    return
end

% ② nargout 已知：按确切个数精确分派
bodyRate = zeros(3, 1);
bodyRateDot = zeros(3, 1);
if numberOfOutputs >= 6
    [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = referenceFcn(t);
    return
end
if numberOfOutputs == 5
    [position, velocity, acceleration, rotation, bodyRate] = referenceFcn(t);
    warnDegraded('5');
    return
end
if numberOfOutputs == 4
    [position, velocity, acceleration, rotation] = referenceFcn(t);
    warnDegraded('4');
    return
end
if numberOfOutputs == 3
    [position, velocity, acceleration] = referenceFcn(t);
    rotation = eye(3);
    warnDegraded('3');
    return
end

% ③ 其余（0~2 个输出）：单输出兜底
[position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    callReferenceSingleOutput(referenceFcn, t, []);
end

function [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    callReferenceAnonymous(referenceFcn, t)
% 匿名句柄专用：nargout == -1，无法预知输出个数，只能逐级回退试探
% 6 → 5 → 4 → 3 个输出，全失败则按单输出解释。
%
% 记下**第一个**错误：若连单输出也失败，就把原始报错抛出去，
% 免得把参考函数内部的真实 bug（维度不符等）伪装成"只返回位置"而静默降级。
bodyRate = zeros(3, 1);
bodyRateDot = zeros(3, 1);
firstError = [];
try
    [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = referenceFcn(t);
    return
catch err
    firstError = err;
end
try
    [position, velocity, acceleration, rotation, bodyRate] = referenceFcn(t);
    warnDegraded('5');
    return
catch
end
try
    [position, velocity, acceleration, rotation] = referenceFcn(t);
    warnDegraded('4');
    return
catch
end
try
    [position, velocity, acceleration] = referenceFcn(t);
    rotation = eye(3);
    warnDegraded('3');
    return
catch
end
[position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    callReferenceSingleOutput(referenceFcn, t, firstError);
end

function [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    callReferenceSingleOutput(referenceFcn, t, firstError)
% 单输出兜底：struct 则读字段，数值则当作位置向量。
bodyRate = zeros(3, 1);
bodyRateDot = zeros(3, 1);
data = [];
try
    data = referenceFcn(t);
catch
    if ~isempty(firstError)
        rethrow(firstError);      % 抛**原始**报错，保留真实病因
    end
    error('crazyflie_slung:referenceFcnUncallable', ...
        'referenceFcn 无法以任何已知契约调用。');
end
if isstruct(data)
    position = data.position;
    velocity = data.velocity;
    acceleration = data.acceleration;
    rotation = data.rotation;
    bodyRate = getFieldOr(data, 'bodyRate', zeros(3, 1));
    bodyRateDot = getFieldOr(data, 'bodyRateDot', zeros(3, 1));
    if ~isfield(data, 'bodyRate') || ~isfield(data, 'bodyRateDot')
        warnDegraded('struct');
    end
else
    position = data;
    velocity = zeros(3, 1);
    acceleration = zeros(3, 1);
    rotation = eye(3);
    warnDegraded('1');
end
end

function warnDegraded(kind)
% 参考函数只提供了部分输出时，一次性提示"哪些量被默认值顶替了"。
% 只警告一次，否则 12500 个积分步会刷满命令窗口。
persistent firedKinds
if isempty(firedKinds)
    firedKinds = {};
end
if any(strcmp(firedKinds, kind))
    return
end
firedKinds{end + 1} = kind;
warning('crazyflie_slung:referenceFcnDegraded', ...
    ['referenceFcn 只提供了 %s 个输出的契约，未提供的量已用默认值补齐' ...
     '（姿态=单位阵 / 体角速度与角加速度=0）。' ...
     '若这不是本意，请检查轨迹函数的输出个数。'], kind);
end

function value = getFieldOr(s, name, fallback)
if isfield(s, name)
    value = s.(name);
else
    value = fallback;
end
end

% ======================================================================
% 日志辅助
% ======================================================================
function c = emptyCommandLog(n)
c = struct(...
    'positionError', zeros(3, 1), ...
    'desiredForce', zeros(3, 1), ...
    'desiredMoment', zeros(3, 1), ...
    'desiredTensions', zeros(3, n), ...
    'tensions', zeros(1, n), ...
    'thrustPercentage', zeros(1, n), ...
    'thrustDesired', zeros(1, n), ...
    'totalThrust', 0, ...
    'desiredLinkUnits', zeros(3, n), ...
    'linkDirectionErrors', zeros(3, n), ...
    'attitudeErrors', zeros(3, n), ...
    'omegaCommands', zeros(3, n), ...
    'omegaCommandDeg', zeros(3, n), ...
    'computedRotations', zeros(3, 3, n));
end

function c = logCommand(command, n)
c = emptyCommandLog(n);
c.positionError = command.positionError;
c.desiredForce = command.desiredForce;
c.desiredMoment = command.desiredMoment;
c.desiredTensions = command.desiredTensions;
c.tensions = command.tensions;
c.thrustPercentage = command.thrustPercentage;
c.thrustDesired = command.thrustDesired;
c.totalThrust = command.totalThrust;
c.desiredLinkUnits = command.desiredLinkUnits;
c.linkDirectionErrors = command.linkDirectionErrors;
c.attitudeErrors = command.attitudeErrors;
c.omegaCommands = command.omegaCommands;
c.omegaCommandDeg = command.omegaCommandDeg;
c.computedRotations = command.computedRotations;
end

% ======================================================================
function summary = computeSummary(sim, cfg)
% 计算关键性能指标，便于自动验证仿真结果是否合理。
nSteps = numel(sim.time);
n = cfg.vehicle.count;
steadyIndex = max(1, round(0.8 * nSteps)):nSteps;

loadPosition = sim.loadPositionLog;
loadHeight = -loadPosition(3, :);
vehicleHeight = -squeeze(sim.vehiclePositionLog(3, :, :));   % n x nSteps

summary.loadPositionFinal = loadPosition(:, end);
summary.loadHeightFinal = loadHeight(end);
summary.vehicleHeightFinal = vehicleHeight(:, end);
summary.steadyPositionError = mean(sim.positionErrorLog(steadyIndex));
summary.maxPositionError = max(sim.positionErrorLog);
summary.steadyLinkError = mean(max(sim.linkErrorLog(:, steadyIndex), [], 1));
summary.maxLinkError = max(sim.linkErrorLog(:));
summary.steadyAttitudeError = mean(max(sim.attitudeErrorLog(:, steadyIndex), [], 1));
summary.maxAttitudeError = max(sim.attitudeErrorLog(:));
% ★ 姿态误差的按轴分解（稳态均值），供自检脚本区分"roll/pitch 可控轴"
%   与"偏航欠驱动轴"。loadAttitudeErrorLog 记录的是每步的三轴姿态误差向量，
%   取稳态窗口的时间均值得到逐轴残差。
%   注意：该日志里各架机的负载姿态误差是同一个量（负载只有一个），
%   因此直接抽出第 1 架即可。
if isfield(sim, 'loadAttitudeErrorLog') && ~isempty(sim.loadAttitudeErrorLog)
    summary.steadyAttitudeErrorVec = mean(sim.loadAttitudeErrorLog(:, steadyIndex), 2);
else
    % 退路：没有按轴日志时用范数标量填充，保证字段存在且维度为 3x1
    summary.steadyAttitudeErrorVec = repmat(summary.steadyAttitudeError, 3, 1);
end
% 负载终态角速度（用于检查偏航漂移率；yaw 分量为欠驱动自由轴）
if isfield(sim, 'loadBodyRateLog') && ~isempty(sim.loadBodyRateLog)
    summary.loadBodyRateFinal = sim.loadBodyRateLog(:, end);
else
    summary.loadBodyRateFinal = zeros(3, 1);
end
summary.maxBodyRate = max(abs(sim.bodyRateLog(:)));
summary.maxThrustPercentage = max(sim.thrustPctLog(:));
summary.minTension = min(sim.tensionLog(:));
summary.maxTension = max(sim.tensionLog(:));
summary.steadyTension = mean(sim.tensionLog(:, steadyIndex), 2);   % n x 1
summary.allTensionsPositive = all(sim.tensionLog(:) > 0);

% -------- 绳索特有断言：长度不变量 ‖挂点 - 无人机‖ ≡ l --------
% ★ 这是"绳（cable）"区别于"刚性连杆（rigid link）"的唯一数值可验量。
%   绳在绷紧时与刚性连杆力学完全相同，长度必须严格守恒；
%   一旦这个量开始漂移，说明 q_i 的约束被破坏、结果不再物理。
%   本实现里无人机位置由 veh = x0 + R0 rho_i - l q_i 定义，
%   因此该不变量**构造性成立**，实测偏差 1.7e-16 m（浮点精度）。
%   这里显式算出来并写进 summary，让自检能守住它，防止将来改动破坏。
maxRopeLengthDrift = 0;
rhoAll = cfg.payload.attachPoints;      % computeSummary 的作用域里没有外部 rhoAll
for k = 1:nSteps
    for i = 1:n
        attachPoint = sim.loadPositionLog(:, k) ...
            + sim.loadRotationLog(:, :, k) * rhoAll(:, i);
        ropeLength = norm(attachPoint - sim.vehiclePositionLog(:, i, k));
        maxRopeLengthDrift = max(maxRopeLengthDrift, ...
            abs(ropeLength - cfg.link.length));
    end
end
summary.maxRopeLengthDrift = maxRopeLengthDrift;
summary.ropeLengthInvariantHolds = maxRopeLengthDrift < 1e-9;

% 每架四旋翼都必须始终高于负载
summary.allVehiclesAboveLoad = all(all(vehicleHeight > loadHeight + 0.05, 2));
summary.finiteState = all(isfinite(loadPosition(:))) ...
    && all(isfinite(sim.linkUnitLog(:))) && all(isfinite(sim.bodyRateLog(:)));
% ★ 坏步计数：> 0 说明代数系统曾经出现非有限量（即使状态日志看起来正常）。
%   自检必须同时要求这一项为 0，否则"发散"会被坏步保护掩盖成"通过"。
if isfield(sim, 'nonFiniteSolveCount')
    summary.nonFiniteSolveCount = sim.nonFiniteSolveCount;
else
    summary.nonFiniteSolveCount = 0;
end
% ★ 发散起点步号（nan = 未触发起点探测）
if isfield(sim, 'onsetStep')
    summary.onsetStep = sim.onsetStep;
else
    summary.onsetStep = nan;
end

% -------- 绕八字避障工况的专项指标 --------
% ★ 三阶段的窗口划分与 cfg.figureEight 完全一致。
%   环绕段（cruise）是"跟踪精度"真正该考核的窗口：起飞/降落段本身就带
%   大范围机动，把它们混进来会把指标稀释得看不出问题。
%   所以把误差统计拆成 起飞 / 环绕 / 降落 三段分别记录。
if isfield(cfg, 'figureEight') && ~isempty(cfg.figureEight)
    fe = cfg.figureEight;
    t1 = fe.takeoffDuration;
    t2 = t1 + fe.cruiseDuration;
    t3 = t2 + fe.landingDuration;
    idxTakeoff = [1, max(2, round(t1 / cfg.simulation.dt))];
    idxCruise = [max(1, round(t1 / cfg.simulation.dt)), ...
                 min(nSteps, round(t2 / cfg.simulation.dt))];
    idxLanding = [min(nSteps, round(t2 / cfg.simulation.dt)), nSteps];
    pe = sim.positionErrorLog;

    summary.phaseWindows = [idxTakeoff; idxCruise; idxLanding];
    summary.maxPositionErrorTakeoff = max(pe(idxTakeoff(1):idxTakeoff(2)));
    summary.maxPositionErrorCruise = max(pe(idxCruise(1):idxCruise(2)));
    summary.meanPositionErrorCruise = mean(pe(idxCruise(1):idxCruise(2)));
    summary.maxPositionErrorLanding = max(pe(idxLanding(1):idxLanding(2)));
    % 降落完成后负载到落点的残差（真正的"能不能停住"指标）
    summary.landingPointError = norm(loadPosition(:, end) - fe.landPosition);
    % 环绕段是"八字"而不是别的路径的判据：x 方向必须完成 2*cycles 次往复
    summary.figureEightCycles = fe.cycles;
    summary.trajectorySpan = [max(loadPosition(1, :)) - min(loadPosition(1, :)); ...
                              max(loadPosition(2, :)) - min(loadPosition(2, :))];
    summary.durationPhases = [t1, t2, t3];
else
    summary.phaseWindows = [];
end

% -------- 锥形障碍物间隙自检 --------
% ★ 论文 Fig. 3 明确要求"around two obstacles represented by cones"，
%   因此"负载是否始终避开锥"是本工况的正确性条件，必须量化守住。
%   判据：负载外接球心到锥轴的**水平**距离，减去锥在该高度处的半径，
%         再减去负载外接球半径，得到净间隙；要求 > cfg.obstacles.clearance。
%   为什么用水平距离而不是三维点面距离：锥是竖直放置的回转体，
%   它的外表面完全由"到轴线的水平距离 r(z)"描述，这样算最直接。
%
% ★★ 高度定义（这里曾经出过一个把整个自检废掉的错误）★★
%   统一用"离地高度"这一个量，避免 baseZ 的符号反复绕：
%       hGround = -loadPosition(3,:)          % z 向下为正，取负即离地高度
%   锥体占据的高度区间是 [hBase, hBase + height]，其中
%       hBase = -cfg.obstacles.baseZ          % 锥底离地高度
%   负载高度落在区间内才比较；高于锥顶则不可能碰撞（gap = inf）。
%   ★ 旧代码写成 hAboveBase = (baseZ - z) + height，再与 height 比，
%     语义混乱：baseZ 为负、height 为正时，实际量是 "离锥顶的距离" 而不是
%     "离锥底的高度"。当锥长期位于负载平面以下时它恒 > height，
%     于是整条轨迹都走 inf 分支、自检恒 PASS。现在改成单一物理量。
hBase = -cfg.obstacles.baseZ;          % 锥底离地高度 [m]
hGround = -loadPosition(3, :);         % 负载离地高度 [m]
if isfield(cfg, 'obstacles') && isfield(cfg.obstacles, 'enabled') ...
        && cfg.obstacles.enabled
    payloadRadius = 0.5 * norm(cfg.payload.size);   % 外接球半径
    nObs = size(cfg.obstacles.positions, 2);
    minClearance = inf;
    minClearanceTime = 0;
    % 诊断：统计各分支被走的次数，防止"恒走 inf 分支"这类静默失效
    nObsInside = 0;
    nObsAbove = 0;
    nObsBelow = 0;
    for k = 1:nSteps
        for j = 1:nObs
            dx = loadPosition(1, k) - cfg.obstacles.positions(1, j);
            dy = loadPosition(2, k) - cfg.obstacles.positions(2, j);
            rHoriz = hypot(dx, dy);
            hLocal = hGround(k) - hBase;          % 相对锥底的高度
            % 只在锥的高度范围内检查；越顶或未进入锥区则间隙视为很大
            if hLocal >= 0 && hLocal <= cfg.obstacles.height
                radiusHere = cfg.obstacles.radius ...
                    * (1 - hLocal / cfg.obstacles.height);
                gap = rHoriz - radiusHere - payloadRadius;
                nObsInside = nObsInside + 1;
            elseif hLocal > cfg.obstacles.height
                gap = inf;    % 负载在锥顶以上，不可能碰撞
                nObsAbove = nObsAbove + 1;
            else
                gap = rHoriz - cfg.obstacles.radius - payloadRadius;
                nObsBelow = nObsBelow + 1;
            end
            if gap < minClearance
                minClearance = gap;
                minClearanceTime = sim.time(k);
            end
        end
    end
    summary.minObstacleClearance = minClearance;
    summary.minObstacleClearanceTime = minClearanceTime;
    summary.obstacleBranchCounts = [nObsInside, nObsAbove, nObsBelow];
    % ★ 有效性门控：若"在锥高度区间内"的采样点数为 0，说明避障根本没被检查过，
    %   此时即便 minClearance 是有限值也毫无意义，必须判定为无效（false）。
    summary.obstacleClearanceValid = (nObsInside > 0);
    summary.obstacleClearanceOk = summary.obstacleClearanceValid ...
        && (minClearance > cfg.obstacles.clearance);
    summary.payloadBoundingRadius = payloadRadius;
else
    summary.minObstacleClearance = inf;
    summary.minObstacleClearanceTime = 0;
    summary.obstacleBranchCounts = [0, 0, 0];
    summary.obstacleClearanceValid = true;
    summary.obstacleClearanceOk = true;
    summary.payloadBoundingRadius = 0.5 * norm(cfg.payload.size);
end

% 悬停时每根绳索张力 = m0 g / n（论文 (18)：sum mu_i = m0 g）
summary.hoverTensionPerLink = cfg.payload.mass * cfg.vehicle.gravity / n;
% 每机悬停推力 = m0 g / n + m_i g（论文 (17) 的悬停解）
summary.hoverThrustPerVehicle = summary.hoverTensionPerLink ...
    + cfg.vehicle.mass * cfg.vehicle.gravity;
summary.hoverThrustTotal = n * summary.hoverThrustPerVehicle;
summary.hoverThrustPercentage = 100 * summary.hoverThrustPerVehicle ...
    / cfg.vehicle.maxTotalThrust;
summary.thrustToWeightRatio = n * cfg.vehicle.maxTotalThrust ...
    / ((cfg.payload.mass + n * cfg.vehicle.mass) * cfg.vehicle.gravity);
end

% ======================================================================
function R = projectSO3(R)
[U, ~, V] = svd(R);
R = U * V.';
if det(R) < 0
    U(:, 3) = -U(:, 3);
    R = U * V.';
end
end

function R = expSO3(v)
theta = norm(v);
if theta < 1e-10
    R = eye(3) + hat(v);
else
    K = hat(v / theta);
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K * K;
end
end

function q = normalizeVector(q)
q = q(:);
n = norm(q);
if n < eps
    q = [0; 0; 1];      % 退化时取 +e3：四旋翼在负载上方（悬停构型）
else
    q = q / n;
end
end

function S = hat(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function value = clamp(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end

function value = clampVector(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end
