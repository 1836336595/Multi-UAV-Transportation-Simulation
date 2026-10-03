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
% 本文件只做：初始化、处理地面起降混合状态、调用控制器、模拟内外环执行机构、
% 积分动力学和记录数据。松弛阶段的无人机位置独立积分；绷紧阶段才使用论文
% 的固定绳长约束。

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
% 初始负载必须满足地面不可穿透约束。z 轴向下为正，所以允许的最大
% z 坐标是 groundZ 减去当前姿态下长方体的竖直半包络高度。
if cfg.takeoff.enabled
    initialGroundLoadZ = cfg.takeoff.groundZ ...
        - payloadGroundHalfHeight(state.loadRotation, cfg.payload.size);
    if state.loadPosition(3) > initialGroundLoadZ
        state.loadPosition(3) = initialGroundLoadZ;
        state.loadVelocity(3) = 0;
    end
end
state.linkUnits = zeros(3, n);
state.linkRates = zeros(3, n);
for i = 1:n
    state.linkUnits(:, i) = normalizeVector(cfg.initial.linkUnits(:, i));
    state.linkRates(:, i) = cfg.initial.linkRates(:, i);
end
state.rotations = zeros(3, 3, n);
state.bodyRates = zeros(3, n);
state.bodyRateDots = zeros(3, n);   % 速率环给出的机体角加速度 Ω̇_i（见下方说明）
for i = 1:n
    state.rotations(:, :, i) = projectSO3(cfg.initial.vehicleR(:, :, i));
    state.bodyRates(:, i) = cfg.initial.bodyRates(:, i);
end
% 松弛绳阶段必须把无人机位置作为独立状态积分；绷紧后再回到论文的
% x_i = x_0 + R_0 rho_i - l_i q_i 约束。这样动画中的绳长不会在起飞时
% 被错误地“瞬间拉满”。
if cfg.takeoff.enabled
    state.vehiclePosition = initialVehiclePositions(state.loadPosition, ...
        state.loadRotation, cfg, false);
else
    state.vehiclePosition = constrainedVehiclePositions(state, cfg);
end
state.vehicleVelocity = zeros(3, n);
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
if cfg.takeoff.enabled
    % 真机起飞阶段电机从零推力开始，由独立位置控制器逐步建立悬停推力；
    % 不能把绷紧悬停解直接作为地面初始推力，否则动画会出现瞬时跳起。
    actualThrust = zeros(1, n);
end

% 控制器跨步状态（与"悬停平衡绳向探针"共用同一份初始化，见 emptyControllerMemory）
memory = emptyControllerMemory(n, cfg);

% 混合状态：独立起飞/收紧 -> 绷紧段 -> 独立释放/降落。
if cfg.takeoff.enabled
    mode = 'SLACK';
else
    mode = 'ACTIVE';
end
modeStartTime = 0;
tautCandidateStart = nan;
tautStartTime = nan;
% ★ 交接瞬间的负载状态快照。TAUT_RAMP 的"参考抬升"必须从这个**冻结**的起点出发，
%   不能再跟着实测跑（见循环内参考抬升一节）。
tautStartPosition = state.loadPosition;
tautStartVelocity = state.loadVelocity;
landingLagWarned = false;
landingStartPosition = state.loadPosition;
groundLoadPosition = state.loadPosition;
groundLoadPosition(3) = cfg.takeoff.groundZ ...
    - payloadGroundHalfHeight(state.loadRotation, cfg.payload.size);
if cfg.takeoff.landingEnabled && cfg.simulation.duration > ...
        cfg.takeoff.takeoffDuration + cfg.takeoff.takeupDuration ...
        + cfg.takeoff.tensionRampTime + cfg.takeoff.landingDuration
    landingStartTime = cfg.simulation.duration - cfg.takeoff.landingDuration;
else
    % 快速冒烟仿真太短时不强行插入降落段，至少保留独立起飞/收紧过程。
    landingStartTime = inf;
end

% ★★★ 收紧段绳向 = **控制器自身的悬停平衡解**（不再是一个统一的几何偏移）
%
% 为什么必须这样（用转储数据定位出来的，不是猜）：
%   `vehicleTakeupPosition` 原来用 `takeoff.groundRadialOffset`（**一个统一值**）
%   把三机放到挂点外侧 ⇒ 三根绳的倾角**完全相同**（本参数下 16.4°）。
%   但悬停张力是 2:1:1，承载大的那根绳更竖直 ⇒ 控制器期望的平衡倾角是三根
%   **不同**的：实测/复算均为 [10.93, 15.27, 15.27]°。
%   于是交接瞬间 `eqi` 一开始就有 5.6/1.2/1.0°，绳向环要在 0.2 s 内吞掉这个阶跃
%   ⇒ 绳向误差冲到 36°、负载角速度 2.9 rad/s、机体速率指令打到限幅。
%   用平衡绳向布置收紧段后，交接瞬间**实际绳向 == 期望绳向 ⇒ eqi ≡ 0**。
%
% ★ 顺带纠正一处概念混淆：`groundRadialOffset` 是**避碰**量
%   （boundingRadius + collisionRadius + clearance），把它当**动力学**的绳向用
%   本来就不成立 —— 两个需求互不相关。地面段照旧用它，收紧段改用平衡绳向。
cfg.link.takeupLinkUnitsBody = equilibriumLinkUnitsBody(state, cfg);

% ★ 安全检查：平衡绳向可能比"地面期的统一外张"更贴近负载（承载大的那根绳更竖直）。
%   这里只用**体系内**几何（与负载朝向无关）算一遍机-负载最小水平距离。
%   ★ 判据用**真实碰撞下限** boundingRadius + collisionRadius，
%     而**不是** 再加上 link.vehicleClearance —— 后者是**地面期**的额外余量，
%     那时机与负载同高、只能靠水平距离保护；收紧时机在负载上方 0.6 m，
%     再把地面余量套上来就又是一次"把两种场景的需求混为一谈"。
%   真低于碰撞下限才告警（不静默放过）。
takeupQBody = cfg.link.takeupLinkUnitsBody;
if ~isempty(takeupQBody)
    takeupDistance = cfg.link.length - cfg.takeoff.preTensionSlack;
    minTakeupHoriz = inf;
    for i = 1:n
        offsetXY = cfg.payload.attachPoints(1:2, i) ...
            - takeupDistance * takeupQBody(1:2, i);
        minTakeupHoriz = min(minTakeupHoriz, norm(offsetXY));
    end
    collisionDistance = cfg.payload.boundingRadius + cfg.vehicle.collisionRadius;
    if minTakeupHoriz < collisionDistance
        warning('crazyflie_slung_simulation:TakeupClearanceTight', ...
            ['收紧段按平衡绳向布置后，机-负载最小水平距离 %.4f m < ' ...
             'boundingRadius + collisionRadius = %.4f m（差 %.4f m）——存在碰撞风险。' ...
             '请调大 allocation.outwardBiasFraction（会同时改变空中平衡绳向），' ...
             '或缩小负载/增大绳长。'], ...
            minTakeupHoriz, collisionDistance, collisionDistance - minTakeupHoriz);
    elseif minTakeupHoriz < collisionDistance + cfg.link.vehicleClearance
        % 只影响地面期的那份额外余量，报告一次（不阻断）
        fprintf(['[takeoff] 提示：收紧段最小机-负载水平距离 %.4f m 已小于' ...
            ' boundingRadius+collisionRadius+vehicleClearance = %.4f m，' ...
            '但仍高于碰撞下限 %.4f m（按平衡绳向布置是有意为之：交接无冲击）。\n'], ...
            minTakeupHoriz, collisionDistance + cfg.link.vehicleClearance, ...
            collisionDistance);
    end
end

% ------------------------------------------------------------------ 日志分配
sim = struct();
sim.time = time;
sim.config = cfg;
sim.loadPositionLog = zeros(3, nSteps);
sim.loadVelocityLog = zeros(3, nSteps);   % ★ 状态量，必须记录（见下面的说明）
sim.loadRotationLog = zeros(3, 3, nSteps);
sim.linkUnitLog = zeros(3, n, nSteps);
sim.vehiclePositionLog = zeros(3, n, nSteps);
sim.ropeDistanceLog = zeros(n, nSteps);
sim.ropeSlackLog = zeros(n, nSteps);
sim.takeoffModeLog = zeros(1, nSteps);
sim.tensionScaleLog = zeros(1, nSteps);
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
sim.bodyRateDotLog = zeros(3, n, nSteps);
sim.attitudeErrorLog = zeros(n, nSteps);
sim.linkErrorLog = zeros(n, nSteps);
sim.positionErrorLog = zeros(1, nSteps);
sim.positionErrorVectorLog = zeros(3, nSteps);
sim.omegaCommandLog = zeros(3, n, nSteps);
sim.desiredTensionLog = zeros(3, n, nSteps);
% 负载姿态误差的按轴日志（3 x nSteps）：用于区分 roll/pitch 与 yaw。
% 负载只有一个，故不重复 n 份。
sim.loadAttitudeErrorLog = zeros(3, nSteps);
% 负载角速度日志（3 x nSteps）：用于检查偏航漂移率。
sim.loadBodyRateLog = zeros(3, nSteps);
% ★★★ 期望负载 yaw 日志（1 x nSteps, rad）= 参考姿态 R0d 第一轴方位角。
%   ★ 定高工况下 R0d ≡ cfg.target.R0（第一轴 = +x）⇒ 本日志**恒为 0**，
%     这是正常数据、不是"没记录到"；画图时**不能**用 any(日志 ~= 0) 判有效
%     （会把恒为 0 的合法数据误判成无数据，于是参考曲线不画）。
sim.loadYawRefLog = zeros(1, nSteps);
% 实际负载 yaw 与包角后的 yaw 跟踪误差
sim.loadYawLog = zeros(1, nSteps);
sim.loadYawErrorLog = zeros(1, nSteps);

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

    % --------------------------- 独立起飞/收紧/释放阶段 ------------------
    % 松弛绳阶段不调用绷紧段张力分配。无人机用 v2 的几何 PID 独立飞行，
    % 负载由地面接触约束固定；绳向 q_i 仅由实际定位几何估计，不需要拉力传感器。
    if cfg.takeoff.enabled && strcmp(mode, 'ACTIVE') ...
            && cfg.takeoff.landingEnabled && t >= landingStartTime
        mode = 'LANDING_TAUT';
        modeStartTime = t;
        landingStartPosition = state.loadPosition;
    elseif cfg.takeoff.enabled && cfg.takeoff.landingEnabled ...
            && t >= landingStartTime ...
            && (strcmp(mode, 'SLACK') || strcmp(mode, 'TAKEUP') ...
            || strcmp(mode, 'TAUT_RAMP'))
        % 若绷紧条件在规定时间内没有满足，也不能让仿真停在“半空等待”。
        % 没有拉力传感器时按保守策略释放并回到地面，日志会保留未绷紧事实。
        mode = 'LANDING_RELEASE';
        modeStartTime = t;
        state.loadPosition = groundLoadPosition;
        state.loadVelocity = zeros(3, 1);
        state.loadBodyRate = zeros(3, 1);
        state.vehicleVelocity = zeros(3, n);
    end

    if cfg.takeoff.enabled && strcmp(mode, 'SLACK') ...
            && t >= cfg.takeoff.takeoffDuration
        mode = 'TAKEUP';
        modeStartTime = t;
    end

    if cfg.takeoff.enabled && (strcmp(mode, 'SLACK') ...
            || strcmp(mode, 'TAKEUP') || strcmp(mode, 'LANDING_RELEASE'))
        if strcmp(mode, 'LANDING_RELEASE')
            state.loadPosition = groundLoadPosition;
            state.loadVelocity = zeros(3, 1);
            state.loadBodyRate = zeros(3, 1);
            [desiredVehicle, desiredVehicleVelocity, desiredVehicleAcceleration] = ...
                independentVehicleTargets(groundLoadPosition, ...
                state.loadRotation, cfg, 'LANDING_RELEASE', t, modeStartTime);
        else
            [desiredVehicle, desiredVehicleVelocity, desiredVehicleAcceleration] = ...
                independentVehicleTargets(groundLoadPosition, ...
                state.loadRotation, cfg, mode, t, modeStartTime);
        end
        [command, memory] = crazyflie_slung_independent_controller(...
            state, desiredVehicle, desiredVehicleVelocity, ...
            desiredVehicleAcceleration, memory, cfg);

        % 独立阶段日志：张力为零，ropeDistance/slack 来自真实位置。
        sim.commandLog(k) = logCommand(command, n);
        sim.tensionLog(:, k) = zeros(n, 1);
        sim.linkErrorLog(:, k) = zeros(n, 1);
        sim.attitudeErrorLog(:, k) = sqrt(sum(command.attitudeErrors.^2, 1)).';
        sim.loadAttitudeErrorLog(:, k) = zeros(3, 1);
        sim.loadBodyRateLog(:, k) = state.loadBodyRate;
        sim.loadYawRefLog(k) = atan2(desired.rotation(2, 1), desired.rotation(1, 1));
        actualYaw = atan2(state.loadRotation(2, 1), state.loadRotation(1, 1));
        sim.loadYawLog(k) = actualYaw;
        sim.loadYawErrorLog(k) = wrapAngle(actualYaw - sim.loadYawRefLog(k));
        sim.positionErrorLog(k) = norm(state.loadPosition - desired.position(:));
        sim.positionErrorVectorLog(:, k) = state.loadPosition - desired.position(:);
        sim.omegaCommandLog(:, :, k) = command.omegaCommands;
        sim.desiredTensionLog(:, :, k) = zeros(3, n);
        sim.loadPositionLog(:, k) = state.loadPosition;
        sim.loadVelocityLog(:, k) = state.loadVelocity;
        sim.loadRotationLog(:, :, k) = state.loadRotation;
        sim.loadBodyRateLog(:, k) = state.loadBodyRate;
        sim.rotationLog(:, :, :, k) = state.rotations;
        sim.bodyRateLog(:, :, k) = state.bodyRates;
        sim.parallelForceLog(:, :, k) = zeros(3, n);
        sim.perpendicularForceLog(:, :, k) = zeros(3, n);
        sim.totalForceLog(:, :, k) = command.totalForces;
        sim.thrustPctLog(:, k) = command.thrustPercentage(:);
        sim.vehiclePositionLog(:, :, k) = state.vehiclePosition;
        sim.takeoffModeLog(k) = modeCode(mode);
        sim.tensionScaleLog(k) = 0;

        if k < nSteps
            % 推力执行器和等效速率环，与绷紧段使用同一套 Crazyflie 模型。
            alphaThrust = min(1, dt / max(cfg.vehicle.thrustTimeConstant, eps));
            actualThrust = actualThrust + alphaThrust * ...
                (command.thrustDesired(:).' - actualThrust);
            actualThrust = min(max(actualThrust, 0), cfg.vehicle.maxTotalThrust);
            sim.thrustLog(:, k + 1) = actualThrust(:);

            % -------- 角速度内环（Crazyflie 固件速率环）：指令 → 角加速度 --------
            % ★★ 与 v2 同构的 **PI** 速率环（不是纯比例）：角速度指令到位后
            %    积分项把稳态速率误差收干净。机体惯量不需要 —— 原来把角加速度包成
            %    力矩 M_i = J_i·α + Ω×J_iΩ 再交给动力学解 Ω̇ = J_i⁻¹(M − Ω×J_iΩ)，
            %    两者恰好抵消 ⇒ 等价于直接给角加速度。
            %    姿态在下面用 expSO3 积分。
            for i = 1:n
                rateError = command.omegaCommands(:, i) - state.bodyRates(:, i);
                memory.rateIntegral(:, i) = memory.rateIntegral(:, i) ...
                    + dt * rateError;
                memory.rateIntegral(:, i) = clampVector( ...
                    memory.rateIntegral(:, i), ...
                    -cfg.rateLoop.integralLimit, cfg.rateLoop.integralLimit);
                state.bodyRateDots(:, i) = ...
                    cfg.rateLoop.bandwidth .* rateError ...
                    + cfg.rateLoop.integralGain .* memory.rateIntegral(:, i);
                sim.bodyRateDotLog(:, i, k + 1) = state.bodyRateDots(:, i);
            end

            uActual = zeros(3, n);
            for i = 1:n
                uActual(:, i) = -actualThrust(i) ...
                    * (state.rotations(:, :, i) * [0; 0; 1]);
                acceleration = cfg.vehicle.gravity * [0; 0; 1] ...
                    + uActual(:, i) / cfg.vehicle.mass;
                state.vehicleVelocity(:, i) = state.vehicleVelocity(:, i) ...
                    + dt * acceleration;
                state.vehiclePosition(:, i) = state.vehiclePosition(:, i) ...
                    + dt * state.vehicleVelocity(:, i);
                % 地面是不可穿透边界。物理坐标 z 向下为正，因此
                % vehicleGroundZ 是机体中心允许达到的最大 z 值。
                vehicleGroundZ = cfg.takeoff.groundZ ...
                    - cfg.takeoff.vehicleGroundClearance;
                if state.vehiclePosition(3, i) > vehicleGroundZ
                    state.vehiclePosition(3, i) = vehicleGroundZ;
                    if state.vehicleVelocity(3, i) > 0
                        state.vehicleVelocity(3, i) = 0;
                    end
                end
                % 机体角速度直接由速率环的角加速度推进（无刚体惯量，见上面的说明）
                state.bodyRates(:, i) = state.bodyRates(:, i) ...
                    + dt * state.bodyRateDots(:, i);
                state.bodyRates(:, i) = clampVector(state.bodyRates(:, i), ...
                    -cfg.simulation.maxBodyRate, cfg.simulation.maxBodyRate);
                state.rotations(:, :, i) = projectSO3(state.rotations(:, :, i) ...
                    * expSO3(state.bodyRates(:, i) * dt));
            end

            [ropeDistance, ropeSlack, qNew, qdNew] = ...
                ropeGeometryFromVehicles(state, cfg);
            state.linkUnits = qNew;
            state.linkRates = qdNew;
            sim.ropeDistanceLog(:, k) = ropeDistance(:);
            sim.ropeSlackLog(:, k) = ropeSlack(:);

            if strcmp(mode, 'TAKEUP')
                % 只有“未超过绳长且余量足够小”才算接近绷直。
                % 仅使用 abs(distance-l) 会把超过绳长的不可实现状态也
                % 当成候选，随后把绳索瞬间投影成刚性约束。
                nearOn = ropeDistance <= cfg.link.length ...
                    & (cfg.link.length - ropeDistance) <= cfg.takeoff.epsilonOn;
                nearOff = ropeDistance <= cfg.link.length ...
                    & (cfg.link.length - ropeDistance) <= cfg.takeoff.epsilonOff;
                if isnan(tautCandidateStart)
                    allNear = all(nearOn);
                else
                    allNear = all(nearOff);
                end
                if allNear
                    if isnan(tautCandidateStart)
                        tautCandidateStart = t;
                    end
                else
                    tautCandidateStart = nan;
                end
                holdTime = t - tautCandidateStart;
                if allNear && holdTime >= cfg.takeoff.confirmTime
                    % ★★ "进入绷紧段时会被模型一次性收掉的绳长余量" = 位置瞬移量。
                    %   本模型只能表示绷紧的绳，所以此刻 |l - d| 会在一拍内被强制归零
                    %   ⇒ 无人机位置瞬移 |l - d|。实测 preTensionSlack = 12 mm 时
                    %   瞬移达 17 mm（等效速率 8.5 m/s），会明显激励绳向环。
                    %   这里显式报出来，不让这种"隐性跳变"悄悄过去。
                    takeupSnap = cfg.link.length - ropeDistance;
                    if max(abs(takeupSnap)) > cfg.takeoff.takeupSnapWarn
                        warning('crazyflie_slung_simulation:TakeupSnap', ...
                            ['进入绷紧段时绳长余量 [%s] m 会被模型一次性收掉，' ...
                             '等效于无人机瞬时位移（最大 %.1f mm）。' ...
                             '把 takeoff.preTensionSlack 取 0 可基本消除' ...
                             '（残余量来自独立控制器的稳态余差）。'], ...
                            mat2str(takeupSnap, 4), ...
                            1000 * max(abs(takeupSnap)));
                    end
                    mode = 'TAUT_RAMP';
                    modeStartTime = t;
                    tautStartTime = t;
                    memory.linkIntegrals = zeros(3, n);
                    memory.previousLinkUnits = state.linkUnits;
                    % 独立段的积分器也清零：交接后再也不用它，留着只会在降落段
                    % 重新启用时带着起飞段的旧偏置（独立控制器已加抗饱和 gate，这里是双保险）。
                    memory.independentPositionIntegral = zeros(3, n);
                    % ★ 冻结交接瞬间的负载状态作为"参考抬升"的起点。
                    tautStartPosition = state.loadPosition;
                    tautStartVelocity = state.loadVelocity;
                end
            end
        end
        continue;
    end

    % ★ 本拍是否处于"交接抬升窗口"。每拍先清零，由下面的抬升块置位；
    %   它同时用于决定**是否允许位置积分器累积**（见控制器调用之后的门控）。
    inHandoffLift = false;

    % LANDING_TAUT 仍使用绷紧段动力学，但把期望负载平滑地送到地面。
    if cfg.takeoff.enabled && strcmp(mode, 'LANDING_TAUT')
        approachDuration = max(0.1, cfg.takeoff.landingDuration ...
            * cfg.takeoff.landingApproachFraction);
        sLanding = clamp((t - modeStartTime) / approachDuration, 0, 1);
        % ★★ 三通道必须同源（同 TAUT_RAMP 的教训）。
        %   旧写法只动了 desired.position，而 desired.velocity / acceleration 恒为 0
        %   ⇒ "参考自己在动、速度参考却是 0"，位置环只能靠反馈硬追 ⇒ 负载明显滞后。
        %   实测（本次转储）：参考 2.2 s 内从 0.350 降到 0.028 m，负载只降到 0.146 m，
        %   滞后 118 mm；随后释放段不得不用"瞬移"把负载按到地面。
        %   改用 5 次多项式剖面：位置/速度/加速度同源，两端速度、加速度均为 0，
        %   与前面的悬停段、后面的释放段都 C^2 连续。
        %   起点用进入降落段时冻结的 landingStartPosition（不再跟实测跑）。
        [blend, blendDot, blendDDot] = smoothStep5WithDerivatives(sLanding);
        dPosLand = groundLoadPosition - landingStartPosition;
        desired.position = landingStartPosition + blend * dPosLand;
        desired.velocity = (blendDot / approachDuration) * dPosLand;
        desired.acceleration = (blendDDot / approachDuration^2) * dPosLand;
        desired.jerk = (smoothStep5ThirdDerivative(sLanding) ...
            / approachDuration^3) * dPosLand;
        desired.rotation = state.loadRotation;
        desired.bodyRate = zeros(3, 1);
        desired.bodyRateDot = zeros(3, 1);
        desired.bodyRateDDot = zeros(3, 1);
    end

    % ★★ 起飞交接：先把张力建起来（参考**冻结**），再把参考按受限剖面抬到目标。
    %
    %   为什么必须分两段（本次转储数据的结论，不再靠猜）：
    %   ① 本设计的悬停总张力**恰好等于负载重量**（hoverTensionByLink 合计 = m0 g），
    %      而 tensionScale 混合的是 uHover = m*g（只抵无人机自重 ⇒ 零张力）。
    %      所以 tensionScale 从 0 到 1 就是张力从 0 到 m0 g
    %      ⇒ **只有斜坡末端（tensionScale ≈ 0.9~1.0）负载才可能离地**。
    %      实测：参考 1.5 s 内从 0.028 升到 0.350 m，而负载到 6.396 s 仍贴地
    %      ⇒ 位置误差在斜坡末端堆到 0.3925 m，之后负载才在 8.7 s 追上。
    %   ② 所以斜坡期间任何"往目标拉"的参考都只会变成纯误差。正确顺序是：
    %      张力建立段（rampT）——参考冻结在交接瞬间的实测状态，误差恒 ≈ 0；
    %      抬升段（liftT）——参考用 5 次多项式剖面从**冻结起点**走到目标，
    %      位置/速度/加速度三通道同源（同下一段的教训），两端速度、加速度均为 0。
    %   ③ 起点必须**冻结**而不能取当前实测值：若锚点是活的，
    %      desired = 实测 + s·(目标 − 实测) ⇒ 误差只能按比例 s 释放，
    %      负载不动时误差照样涨满（这正是上一版的行为）。
    if cfg.takeoff.enabled && ~isnan(tautStartTime)
        rampT = max(cfg.takeoff.tensionRampTime, eps);
        liftT = max(cfg.takeoff.referenceLiftTime, eps);
        tSinceTaut = t - tautStartTime;
        inLiftWindow = (tSinceTaut >= 0) && (tSinceTaut < rampT + liftT) ...
            && (strcmp(mode, 'TAUT_RAMP') || strcmp(mode, 'ACTIVE'));
        if inLiftWindow
            inHandoffLift = true;
            sLift = clamp((tSinceTaut - rampT) / liftT, 0, 1);
            [refBlend, refBlendDot, refBlendDDot] = ...
                smoothStep5WithDerivatives(sLift);
            posTarget = desired.position;
            velTarget = desired.velocity;
            accTarget = desired.acceleration;
            jerkTarget = desired.jerk;
            dPos = posTarget - tautStartPosition;
            dVel = velTarget - tautStartVelocity;
            desired.position     = tautStartPosition + refBlend * dPos;
            desired.velocity     = tautStartVelocity ...
                + (refBlendDot / liftT) * dPos + refBlend * dVel;
            desired.acceleration = (refBlendDDot / liftT^2) * dPos ...
                + refBlend * accTarget;
            desired.jerk = (smoothStep5ThirdDerivative(sLift) / liftT^3) * dPos ...
                + (refBlendDDot / liftT^2) * velTarget ...
                + (refBlendDot / liftT) * accTarget + refBlend * jerkTarget;
        end
    end

    % -------- 控制器 --------
    [command, memory] = crazyflie_slung_controller(state, desired, memory, cfg);

    % ★★ 非绷紧阶段禁止负载位置积分器累积（2026-09-28 修复）
    %   机理：TAUT_RAMP 期间控制器照常被调用（斜坡只缩放它的**输出**，见下方
    %   tensionScale），而此刻负载还贴着地面、离目标高度误差很大
    %   ⇒ positionIntegral 一路冲到限幅 0.5，斜坡结束时把『憋住的力』
    %     一次性释放 ⇒ 交接峰值反而更高。
    %   ★ 实证：把 tensionRampTime 从 1.5 s 拉长到 5.5 s 后，绳向误差峰值
    %     由 26 deg 恶化到 50~60 deg —— 正是这条机理（斜坡越长憋得越久）。
    %   所以只在绷紧运输阶段（ACTIVE / LANDING_TAUT）让它累积。
    %   ★★★ 2026-09-28 修复：这里原写成 'TAUT_ACTIVE'，但状态机里的名字是
    %     'ACTIVE'（见文件末尾 modeCode 的 case）⇒ 条件恒真 ⇒
    %     **整个运输段每一步都把 positionIntegral 清零，ki 完全失效**。
    %     本项目里 'TAUT_ACTIVE' 只是文档/图例里的叫法（README、
    %     visualization 的标签），代码里的实际字符串是 'ACTIVE'。
    %     ⇒ 凡是用 strcmp(mode, ...) 的地方，名字**必须**取自 modeCode 的 case。
    %   ★ 另外，交接抬升窗口（inHandoffLift）里也要清零：那一段是"指令性瞬态"，
    %     参考由我们自己的前馈剖面给出，残余误差不是常值扰动；若让积分器累积，
    %     抬升结束时它会带着一份"憋住的力"把负载顶过目标（新的过冲来源）。
    %   ★ 注意：一旦 'ACTIVE' 的名字改对，ki 就真正生效了 —— 巡航段的行为
    %     会与之前（ki 被无意清零）不同，需要重新确认。
    if ~((strcmp(mode, 'ACTIVE') || strcmp(mode, 'LANDING_TAUT')) && ~inHandoffLift)
        memory.positionIntegral = zeros(3, 1);
    end

    % 绷紧瞬间不要把完整张力阶跃施加到地面负载。把各机期望作用力从
    % 独立悬停力平滑过渡到论文绷紧段作用力；这只是接触阶段的数值/物理
    % 过渡，ACTIVE 段中仍完全使用原控制律。
    tensionScale = 1.0;
    if strcmp(mode, 'TAUT_RAMP')
        rampS = clamp((t - tautStartTime) / ...
            max(cfg.takeoff.tensionRampTime, eps), 0, 1);
        tensionScale = smoothStep5(rampS);
        uHover = -cfg.vehicle.mass * cfg.vehicle.gravity * [0; 0; 1];
        uBlend = repmat(uHover, 1, n) ...
            + tensionScale * (command.totalForces - repmat(uHover, 1, n));
        command.totalForces = uBlend;
        command.parallelForces = tensionScale * command.parallelForces ...
            + (1 - tensionScale) * repmat(uHover, 1, n);
        command.perpendicularForces = tensionScale * command.perpendicularForces;
        command.desiredTensions = tensionScale * command.desiredTensions;
        command.tensions = tensionScale * command.tensions;
        command.thrustDesired = zeros(1, n);
        for i = 1:n
            command.thrustDesired(i) = clamp(-dot(uBlend(:, i), ...
                state.rotations(:, :, i) * [0; 0; 1]), ...
                0, cfg.vehicle.maxTotalThrust);
        end
        command.thrustPercentage = 100 * command.thrustDesired ...
            / cfg.vehicle.maxTotalThrust;
        command.totalThrust = sum(command.thrustDesired);
    end
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
    sim.loadYawRefLog(k) = atan2(desired.rotation(2, 1), desired.rotation(1, 1));
    actualYaw = atan2(state.loadRotation(2, 1), state.loadRotation(1, 1));

    referenceYaw = atan2(desired.rotation(2, 1), desired.rotation(1, 1));

    sim.loadYawLog(k) = actualYaw;
    sim.loadYawErrorLog(k) = wrapAngle(actualYaw - referenceYaw);
    sim.positionErrorLog(k) = norm(command.positionError);
    sim.positionErrorVectorLog(:, k) = command.positionError;
    sim.omegaCommandLog(:, :, k) = command.omegaCommands;
    sim.desiredTensionLog(:, :, k) = command.desiredTensions;
    sim.takeoffModeLog(k) = modeCode(mode);
    sim.tensionScaleLog(k) = tensionScale;

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
        state.vehiclePosition(:, i) = state.loadPosition ...
            + state.loadRotation * rhoAll(:, i) ...
            - cfg.link.length * state.linkUnits(:, i);
        sim.vehiclePositionLog(:, i, k) = state.vehiclePosition(:, i);
        attachPoint = state.loadPosition + state.loadRotation * rhoAll(:, i);
        sim.ropeDistanceLog(i, k) = norm(attachPoint - state.vehiclePosition(:, i));
        sim.ropeSlackLog(i, k) = cfg.link.length - sim.ropeDistanceLog(i, k);
    end

    if k == nSteps
        break;
    end

    % -------- 推力执行器：一阶响应（逐机独立） --------
    alpha = min(1, dt / max(cfg.vehicle.thrustTimeConstant, eps));
    actualThrust = actualThrust + alpha * (command.thrustDesired(:).' - actualThrust);
    actualThrust = min(max(actualThrust, 0), cfg.vehicle.maxTotalThrust);
    sim.thrustLog(:, k + 1) = actualThrust(:);

    % -------- 角速度内环：Crazyflie 固件速率环（逐机独立） --------
    % ★ 与 v2 同构的 **PI**：Omega_dot_i = bandwidth.*e_i + integralGain.*∫e_i，
    %   e_i = Omega_cmd_i - Omega_i。只给角加速度，**不再包成力矩**
    %   （机体惯量会与姿态方程里的叉乘项完全抵消，见 parameters.m 的推导）。
    for i = 1:n
        rateError = command.omegaCommands(:, i) - state.bodyRates(:, i);
        memory.rateIntegral(:, i) = memory.rateIntegral(:, i) ...
            + dt * rateError;
        memory.rateIntegral(:, i) = clampVector(memory.rateIntegral(:, i), ...
            -cfg.rateLoop.integralLimit, cfg.rateLoop.integralLimit);
        state.bodyRateDots(:, i) = cfg.rateLoop.bandwidth .* rateError ...
            + cfg.rateLoop.integralGain .* memory.rateIntegral(:, i);
        sim.bodyRateDotLog(:, i, k + 1) = state.bodyRateDots(:, i);
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

    % 地面接触和混合状态切换。负载中心不能低于“底面接触地面”的位置；
    % 这一步避免绷紧/降落阶段出现数值穿地，也让可视化中的接触过程可信。
    if cfg.takeoff.enabled
        % z 轴向下为正：超过该值表示负载底面已经穿过地面。
        currentGroundLoadZ = cfg.takeoff.groundZ ...
            - payloadGroundHalfHeight(state.loadRotation, cfg.payload.size);
        if state.loadPosition(3) > currentGroundLoadZ
            state.loadPosition(3) = currentGroundLoadZ;
            if state.loadVelocity(3) > 0
                state.loadVelocity(3) = 0;
            end
        end

        % ★★ 地面摩擦（2026-09-29 新增）—— 竖直方向夹住了，**水平方向也必须给摩擦**
        %
        % 原来这里只有 z 向约束、没有 xy 向摩擦 ⇒ 负载"躺在无摩擦地面上"。
        % 实测后果（转储数据）：交接期负载高度**一直贴在 0.0279 m（就是触地高度）**，
        % 却被绳子的水平不平衡力拖走 **0.18 m**（err_x 5.778 s →0，6.9 s →−0.179 m）。
        % 这条"持续扰动"又反过来改变挂点位置、激励绳向环：
        % 绳向误差冲到 30°、负载角速度 2.7 rad/s、速率/姿态指令打到限幅。
        % 但真实情况下静摩擦上限 ~μ*m0*g = 0.24~0.39 N，比那个扰动（~0.06 N）大
        % 3~6 倍 ⇒ **负载根本不会滑**。所以原来的滑动是**模型缺摩擦**造的假象。
        %
        % 实现：库仑摩擦（水平减速度上限 μ*g），且 μ 随**法向力**衰减 ——
        % 绳把负载往上提 ⇒ 法向力减小 ⇒ 摩擦上限减小；
        % 张力涨到等于自重时法向力→0、摩擦自然消失，不会把负载"粘"在半空。
        % 这样离地瞬间是干净的，不需要额外的释放判据。
        if state.loadPosition(3) >= currentGroundLoadZ ...
                - cfg.takeoff.groundContactTolerance
            % ★ 注意本文件里没有 m0 / g 这两个简写（只有 cfg.payload.mass / cfg.vehicle.gravity），
            %   写 m0*g 会直接报"函数或变量无法识别"。
            weight = cfg.payload.mass * cfg.vehicle.gravity;
            normalLoad = weight ...
                - sum(command.tensions(:).' .* state.linkUnits(3, :));
            frictionRatio = max(0, min(1, normalLoad / weight));
            frictionMu = cfg.takeoff.groundFrictionMu * frictionRatio;
            horizontalSpeed = norm(state.loadVelocity(1:2));
            if horizontalSpeed > 0 && frictionMu > 0
                speedDrop = min(frictionMu * cfg.vehicle.gravity * dt, ...
                    horizontalSpeed);
                state.loadVelocity(1:2) = state.loadVelocity(1:2) ...
                    * (1 - speedDrop / horizontalSpeed);
            end

            % ★★ 地面也**约束转动**（2026-09-29 补）—— 这是交接期"绳向误差 30°"的真源头
            %
            % 躺在地面上的负载不可能以 ~3 rad/s 翻滚 —— 那要求把一条边抬起来。
            % 但原来这里只约束 z、完全不约束转动，于是交接期负载在地面上"打滚"：
            % 实测 |ω| 达 2.9 rad/s，分量主要在 ω_x/ω_y（ω_z 只有 ~0.1）。
            %
            % 后果不只是不真实 —— 角速度经
            %       Md = -kR*eR0 - kOmega*eOmega0 + ...
            % 直接进**张力分配**：kOmega = [0.00276, 0.00403, 0.000339]，
            % ω = 2.9 rad/s ⇒ |Md| 可达 0.012 N·m；而力臂只有 ~0.04 m
            % ⇒ 等效张力扰动 0.012/(3*0.04) ≈ 0.10 N —— 是绳 2 张力(0.20 N)的一半！
            % ⇒ **期望绳向 q_id 被甩来甩去**，`eqi = cross(q_id, q_i)` 冲到 30°，
            %   而实际绳向其实只偏 1~2°（用无人机/负载位置重建可验证：
            %   `_verify_tools/_diag_link_azimuth.py` 给出重算夹角 ≤2.3°）。
            % ⇒ 那个"30° 绳向误差"主要是**被地面上的打滚甩出来的假象**。
            %
            % 摩擦转矩上限 ≈ mu*m0*g*rho_typ，除以 J0 得角减速度上限：
            %    0.5 * 0.08 * 9.81 * 0.04 / 2.69e-4 ≈ 58 rad/s²（滚转/俯仰）
            % 取 40 rad/s^2 已足够锁住触地期的转动；同样随法向力衰减，
            % 张力把负载提起来时制动自然消失（与水平摩擦一致）。
            brakeRate = cfg.takeoff.groundRotationalBrake * frictionRatio;
            bodyRateNorm = norm(state.loadBodyRate);
            if bodyRateNorm > 0 && brakeRate > 0
                rateDrop = min(brakeRate * dt, bodyRateNorm);
                state.loadBodyRate = state.loadBodyRate ...
                    * (1 - rateDrop / bodyRateNorm);
            end
        end
    end

    if cfg.takeoff.enabled && strcmp(mode, 'LANDING_TAUT')
        approachDuration = max(0.1, cfg.takeoff.landingDuration ...
            * cfg.takeoff.landingApproachFraction);
        % ★★ 释放条件：时刻到了**并且**负载已经落到地面附近。
        %   旧写法只判时刻 ⇒ 若负载落后于参考（本次实测落后 118 mm），
        %   释放瞬间 state.loadPosition = groundLoadPosition 会把负载**瞬移**到地面，
        %   日志上是高度阶跃、物理上是不存在的卸载冲击，而且**不会报错**。
        %   ⇒ 改成"落到地面附近才释放"，并保证还有绷紧动力学把它压下来；
        %     超时未落地时不瞬移，只发一次 warning（静默失效是本项目最大的坑）。
        %   z 轴向下为正 ⇒ 负载"高于地面"时 z < currentGroundLoadZ，
        %   离地高度 = currentGroundLoadZ - z。
        nearGround = state.loadPosition(3) >= currentGroundLoadZ ...
            - cfg.takeoff.releaseSnapTolerance;
        timeUp = t >= modeStartTime + approachDuration;
        if timeUp && nearGround
            state.loadPosition = groundLoadPosition;
            state.loadVelocity = zeros(3, 1);
            state.loadBodyRate = zeros(3, 1);
            mode = 'LANDING_RELEASE';
            modeStartTime = t;
            state.vehiclePosition = constrainedVehiclePositions(state, cfg);
            state.vehicleVelocity = zeros(3, n);
            [~, ~, state.linkUnits, state.linkRates] = ...
                ropeGeometryFromVehicles(state, cfg);
        elseif timeUp && ~landingLagWarned
            landingLagWarned = true;
            warning('crazyflie_slung_simulation:LandingApproachLag', ...
                ['下降参考已到地面，但负载仍离地 %.1f mm（> releaseSnapTolerance %.1f mm）。' ...
                 '继续用绷紧动力学下降，不把负载瞬移到地面。'], ...
                1000 * (currentGroundLoadZ - state.loadPosition(3)), ...
                1000 * cfg.takeoff.releaseSnapTolerance);
        end
    elseif cfg.takeoff.enabled && strcmp(mode, 'TAUT_RAMP') ...
            && t >= tautStartTime + cfg.takeoff.tensionRampTime
        if state.loadPosition(3) > groundLoadPosition(3)
            state.loadPosition(3) = groundLoadPosition(3);
            state.loadVelocity(3) = 0;
        end
        mode = 'ACTIVE';
        modeStartTime = t;
    elseif cfg.takeoff.enabled && strcmp(mode, 'TAUT_RAMP') ...
            && state.loadPosition(3) > groundLoadPosition(3)
        state.loadPosition(3) = groundLoadPosition(3);
        state.loadVelocity(3) = 0;
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
    if cfg.takeoff.enabled && (strcmp(mode, 'SLACK') ...
            || strcmp(mode, 'TAKEUP') || strcmp(mode, 'LANDING_RELEASE'))
        state.vehiclePosition(:, i) = state.vehiclePosition(:, i);
    else
        state.vehiclePosition(:, i) = state.loadPosition ...
            + state.loadRotation * rhoAll(:, i) ...
            - cfg.link.length * state.linkUnits(:, i);
    end
    sim.vehiclePositionLog(:, i, nSteps) = state.vehiclePosition(:, i);
    attachPoint = state.loadPosition + state.loadRotation * rhoAll(:, i);
    sim.ropeDistanceLog(i, nSteps) = norm(attachPoint - state.vehiclePosition(:, i));
    sim.ropeSlackLog(i, nSteps) = cfg.link.length - sim.ropeDistanceLog(i, nSteps);
end
sim.takeoffModeLog(nSteps) = modeCode(mode);

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
% 读取静态目标。本工程只有**定高**工况：目标是常量，所有导数通道都有确定取值。
%
% ★ 2026-10-03：这里原先还有"`cfg.referenceFcn` 非空则调用外部轨迹函数"的分支，
%   连同 8/7/6/5/4/3 输出的多契约兼容层（`callReferenceFunction` 等）一起删掉了 ——
%   八字轨迹只存在于另一个分支，本分支只做定高。
%   交接抬升 / 下降段的**时变**参考由主循环里的 5 次多项式剖面覆盖（见调用处），
%   它的 jerk 是解析算出来的，所以 analytic 前馈要的那两个量依然齐备。
%   ⇒ 顺带说明：原来那段"analytic 必须有 jerk/bodyRateDDot"的严格检查也随之删除，
%     因为它现在**不可能触发**（静态路径直接给 0，抬升/降落段覆盖成解析值）。
desired.position = cfg.target.position;
desired.velocity = cfg.target.velocity;
desired.acceleration = cfg.target.acceleration;
desired.rotation = cfg.target.R0;
desired.bodyRate = cfg.reference.omegaD;
desired.bodyRateDot = cfg.reference.omegaDotD;
desired.jerk = zeros(3, 1);          % 静态目标 ⇒ 位置高阶导数恒为 0
desired.bodyRateDDot = zeros(3, 1);  % 偏航锁定时负载期望角速度恒为 0，其导数亦然

desired.position = desired.position(:);
desired.velocity = desired.velocity(:);
desired.acceleration = desired.acceleration(:);
desired.rotation = projectSO3(desired.rotation);
desired.bodyRate = desired.bodyRate(:);
desired.bodyRateDot = desired.bodyRateDot(:);
desired.jerk = desired.jerk(:);
desired.bodyRateDDot = desired.bodyRateDDot(:);
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
% 倾斜缆绳不要求无人机的水平投影落在负载内部；这里只记录垂直净空，
% 用于诊断无人机是否仍位于负载平面上方，而不把它误当成“正上方”约束。
if isvector(vehicleHeight)
    vehicleHeight = reshape(vehicleHeight, n, nSteps);
end

summary.loadPositionFinal = loadPosition(:, end);
summary.loadHeightFinal = loadHeight(end);
summary.vehicleHeightFinal = vehicleHeight(:, end);
if isfield(sim, 'takeoffModeLog')
    summary.takeoffModeFinal = sim.takeoffModeLog(end);
    summary.takeoffModeCodesSeen = unique(sim.takeoffModeLog);
    summary.tautTransitionOccurred = any(sim.takeoffModeLog == 2 ...
        | sim.takeoffModeLog == 3 | sim.takeoffModeLog == 4);
    summary.firstTautTime = NaN;
    firstTautIndex = find(sim.takeoffModeLog == 2 | sim.takeoffModeLog == 3 ...
        | sim.takeoffModeLog == 4, 1, 'first');
    if ~isempty(firstTautIndex)
        summary.firstTautTime = sim.time(firstTautIndex);
    end
    summary.finalRopeDistance = sim.ropeDistanceLog(:, end);
    summary.finalRopeSlack = sim.ropeSlackLog(:, end);
else
    summary.takeoffModeFinal = 3;
    summary.takeoffModeCodesSeen = 3;
    summary.tautTransitionOccurred = true;
    summary.firstTautTime = 0;
    summary.finalRopeDistance = repmat(cfg.link.length, n, 1);
    summary.finalRopeSlack = zeros(n, 1);
end
if cfg.takeoff.enabled
    vehicleGroundZ = cfg.takeoff.groundZ - cfg.takeoff.vehicleGroundClearance;
    summary.maxVehicleGroundPenetration = max(0, max(sim.vehiclePositionLog(3, :, :), [], 'all') ...
        - vehicleGroundZ);
    payloadGroundZ = zeros(1, nSteps);
    for k = 1:nSteps
        payloadGroundZ(k) = cfg.takeoff.groundZ ...
            - payloadGroundHalfHeight(sim.loadRotationLog(:, :, k), ...
            cfg.payload.size);
    end
    summary.maxPayloadGroundPenetration = max(0, ...
        max(loadPosition(3, :) - payloadGroundZ));
else
    summary.maxVehicleGroundPenetration = 0;
    summary.maxPayloadGroundPenetration = 0;
end
% yaw 误差必须逐点包角后再统计，避免 ±pi 处产生假大误差
if isfield(sim, 'loadYawErrorLog')
    yawError = sim.loadYawErrorLog(:).';

    summary.steadyYawTrackingError = ...
        mean(abs(yawError(steadyIndex)));

    summary.maxYawTrackingError = ...
        max(abs(yawError));

    summary.finalYawTrackingError = ...
        abs(yawError(end));

    summary.loadYawFinal = sim.loadYawLog(end);
    summary.loadYawReferenceFinal = sim.loadYawRefLog(end);
else
    summary.steadyYawTrackingError = NaN;
    summary.maxYawTrackingError = NaN;
    summary.finalYawTrackingError = NaN;
    summary.loadYawFinal = NaN;
    summary.loadYawReferenceFinal = NaN;
end
summary.steadyPositionError = mean(sim.positionErrorLog(steadyIndex));
summary.maxPositionError = max(sim.positionErrorLog);
summary.steadyLinkError = mean(max(sim.linkErrorLog(:, steadyIndex), [], 1));
summary.maxLinkError = max(sim.linkErrorLog(:));
summary.steadyAttitudeError = mean(max(sim.attitudeErrorLog(:, steadyIndex), [], 1));
summary.maxAttitudeError = max(sim.attitudeErrorLog(:));
% ★ 姿态误差的按轴分解（稳态均值），供自检脚本区分"roll/pitch 可控轴"
%   与 yaw 弱可控轴。loadAttitudeErrorLog 记录的是每步的三轴姿态误差向量，
%   取稳态窗口的时间均值得到逐轴残差。
%   注意：该日志里各架机的负载姿态误差是同一个量（负载只有一个），
%   因此直接抽出第 1 架即可。
if isfield(sim, 'loadAttitudeErrorLog') && ~isempty(sim.loadAttitudeErrorLog)
    summary.steadyAttitudeErrorVec = mean(sim.loadAttitudeErrorLog(:, steadyIndex), 2);
else
    % 退路：没有按轴日志时用范数标量填充，保证字段存在且维度为 3x1
    summary.steadyAttitudeErrorVec = repmat(summary.steadyAttitudeError, 3, 1);
end
% 负载终态角速度（用于检查 yaw 闭环的角速度收敛）
if isfield(sim, 'loadBodyRateLog') && ~isempty(sim.loadBodyRateLog)
    summary.loadBodyRateFinal = sim.loadBodyRateLog(:, end);
else
    summary.loadBodyRateFinal = zeros(3, 1);
end
summary.maxBodyRate = max(abs(sim.bodyRateLog(:)));
summary.maxThrustPercentage = max(sim.thrustPctLog(:));
summary.minTension = min(sim.tensionLog(:));
summary.maxTension = max(sim.tensionLog(:));
tautMask = true(1, nSteps);
if isfield(sim, 'takeoffModeLog')
    tautMask = sim.takeoffModeLog == 2 | sim.takeoffModeLog == 3 ...
        | sim.takeoffModeLog == 4;
end
summary.hasTautPhase = any(tautMask);
positiveTensionMask = sim.takeoffModeLog == 3 | sim.takeoffModeLog == 4;
if any(positiveTensionMask)
    summary.minTautTension = min(sim.tensionLog(:, positiveTensionMask), [], 'all');
    summary.maxTautTension = max(sim.tensionLog(:, positiveTensionMask), [], 'all');
else
    summary.minTautTension = 0;
    summary.maxTautTension = 0;
end
steadyTensionIndex = steadyIndex(tautMask(steadyIndex));
if isempty(steadyTensionIndex)
    steadyTensionIndex = find(tautMask);
end
if isempty(steadyTensionIndex)
    summary.steadyTension = zeros(n, 1);
else
    summary.steadyTension = mean(sim.tensionLog(:, steadyTensionIndex), 2);
end
if any(positiveTensionMask)
    summary.allTensionsPositive = all(sim.tensionLog(:, positiveTensionMask) > 0, 'all');
else
    summary.allTensionsPositive = true;
end

% 无人机碰撞诊断：用机臂长度+旋翼半径作为保守水平包络半径。
% 绳索允许倾斜后，车辆中心不再被强制固定在挂点正上方；这里直接检查
% 记录到的实际车辆中心距，避免“几何上有倾斜、图上仍发生机体重叠”。
minVehicleSeparation = inf;
for k = 1:nSteps
    for i = 1:n
        for j = i + 1:n
            centerDistance = norm(sim.vehiclePositionLog(:, i, k) ...
                - sim.vehiclePositionLog(:, j, k));
            minVehicleSeparation = min(minVehicleSeparation, centerDistance);
        end
    end
end
summary.minVehicleSeparation = minVehicleSeparation;
requiredVehicleSeparation = 2 * cfg.vehicle.collisionRadius ...
    + cfg.link.vehicleClearance;
summary.requiredVehicleSeparation = requiredVehicleSeparation;
summary.vehicleCollisionFree = minVehicleSeparation > requiredVehicleSeparation;
summary.vehicleSeparationMargin = minVehicleSeparation - requiredVehicleSeparation;
% 无人机-负载碰撞诊断：用负载外接球和无人机水平包络半径作保守估计。
% 该间隙随 payload.size 自动变化；它只用于诊断，不改变绳索动力学。
minVehiclePayloadClearance = inf;
for k = 1:nSteps
    for i = 1:n
        centerDistance = norm(sim.vehiclePositionLog(:, i, k) - loadPosition(:, k));
        gap = centerDistance - cfg.payload.boundingRadius ...
            - cfg.vehicle.collisionRadius;
        minVehiclePayloadClearance = min(minVehiclePayloadClearance, gap);
    end
end
summary.minVehiclePayloadClearance = minVehiclePayloadClearance;
summary.vehiclePayloadCollisionFree = ...
    minVehiclePayloadClearance > cfg.link.vehicleClearance;
summary.vehiclePayloadClearanceMargin = minVehiclePayloadClearance ...
    - cfg.link.vehicleClearance;

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
    if ~tautMask(k)
        continue;
    end
    for i = 1:n
        attachPoint = sim.loadPositionLog(:, k) ...
            + sim.loadRotationLog(:, :, k) * rhoAll(:, i);
        ropeLength = norm(attachPoint - sim.vehiclePositionLog(:, i, k));
        maxRopeLengthDrift = max(maxRopeLengthDrift, ...
            abs(ropeLength - cfg.link.length));
    end
end
summary.maxRopeLengthDrift = maxRopeLengthDrift;
summary.ropeLengthInvariantHolds = ~summary.hasTautPhase ...
    || maxRopeLengthDrift < 1e-9;

% 倾斜缆绳下，车辆可以明显偏离负载的水平投影，不能再要求“在负载正上方”。
% 保留旧字段以兼容外部脚本，但它现在只表示垂直净空诊断，不参与控制律。
verticalClearance = vehicleHeight - repmat(loadHeight, n, 1);
summary.minVehicleVerticalClearance = min(verticalClearance(:));
summary.allVehiclesAboveLoad = summary.minVehicleVerticalClearance > 0;
summary.minLinkVerticalComponent = min(sim.linkUnitLog(3, :, :), [], 'all');
summary.finiteState = all(isfinite(loadPosition(:))) ...
    && all(isfinite(sim.linkUnitLog(:))) && all(isfinite(sim.bodyRateLog(:))) ...
    && all(isfinite(sim.vehiclePositionLog(:))) ...
    && all(isfinite(sim.ropeDistanceLog(:)));
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


% 悬停张力由当前挂点几何决定，而不是固定的 m0*g/n。
% 对竖直绳索，sum(mu_i) = -m0*g*e3 且 sum(rho_i x mu_i) = 0，
% 因此用当前 rho_i 解出每根绳的解析悬停张力。尺寸改变后，力臂和该向量
% 会自动更新；这也与 demo 中的分配一致性检查使用同一几何。
rhoHover = cfg.payload.attachPoints;
hoverMatrix = [ones(1, n); rhoHover(2, :); -rhoHover(1, :)];
if size(hoverMatrix, 1) <= size(hoverMatrix, 2) ...
        && rank(hoverMatrix) == size(hoverMatrix, 1)
    signedHoverTension = hoverMatrix \ [-cfg.payload.mass * cfg.vehicle.gravity; 0; 0];
    summary.hoverTensionByLink = abs(signedHoverTension(:));
else
    summary.hoverTensionByLink = repmat(cfg.payload.mass * cfg.vehicle.gravity / n, n, 1);
end
summary.hoverTensionPerLink = mean(summary.hoverTensionByLink);
% 每机悬停推力 = 当前绳张力 + m_i*g（论文 (17) 的竖直悬停解）
summary.hoverThrustByVehicle = summary.hoverTensionByLink ...
    + cfg.vehicle.mass * cfg.vehicle.gravity;
summary.hoverThrustPerVehicle = mean(summary.hoverThrustByVehicle);
summary.hoverThrustTotal = sum(summary.hoverThrustByVehicle);
summary.hoverThrustPercentage = 100 * summary.hoverThrustByVehicle ...
    / cfg.vehicle.maxTotalThrust;
summary.thrustToWeightRatio = n * cfg.vehicle.maxTotalThrust ...
    / ((cfg.payload.mass + n * cfg.vehicle.mass) * cfg.vehicle.gravity);
end

% ======================================================================
function positions = initialVehiclePositions(loadPosition, loadRotation, cfg, useTakeupHeight)
% 初始时无人机在地面附近，绳索必然松弛；随后由独立控制器起飞。
n = cfg.vehicle.count;
rhoAll = cfg.payload.attachPoints;
positions = zeros(3, n);
for i = 1:n
    attach = loadPosition + loadRotation * rhoAll(:, i);
    radial = horizontalRadial(rhoAll(:, i), i, n);
    positions(:, i) = attach + cfg.takeoff.groundRadialOffset * radial;
    if useTakeupHeight
        positions(:, i) = vehicleTakeupPosition(attach, radial, loadRotation, i, cfg);
    else
        positions(3, i) = cfg.takeoff.groundZ ...
            - cfg.takeoff.vehicleGroundClearance;
    end
end
end

function positions = constrainedVehiclePositions(state, cfg)
n = cfg.vehicle.count;
positions = zeros(3, n);
for i = 1:n
    positions(:, i) = state.loadPosition ...
        + state.loadRotation * cfg.payload.attachPoints(:, i) ...
        - cfg.link.length * state.linkUnits(:, i);
end
end

function [positions, velocities, accelerations] = ...
    independentVehicleTargets(loadPosition, loadRotation, cfg, mode, t, modeStartTime)
% 返回每架独立飞行阶段的目标位置。TAKEUP 目标是 l-preTensionSlack，
% 因此不会在尚未确认所有绳索接近绷直时提前切换到耦合动力学。
n = cfg.vehicle.count;
rhoAll = cfg.payload.attachPoints;
positions = zeros(3, n);
velocities = zeros(3, n);
accelerations = zeros(3, n);
for i = 1:n
    attach = loadPosition + loadRotation * rhoAll(:, i);
    radial = horizontalRadial(rhoAll(:, i), i, n);
    if strcmp(mode, 'TAKEUP')
        takeupPosition = vehicleTakeupPosition(attach, radial, loadRotation, i, cfg);
        hoverPosition = attach + cfg.takeoff.groundRadialOffset * radial;
        hoverPosition(3) = cfg.takeoff.groundZ ...
            - cfg.takeoff.independentHoverHeight;
        takeupDuration = max(cfg.takeoff.takeupDuration, eps);
        takeupU = (t - modeStartTime) / takeupDuration;
        [takeupBlend, takeupBlendDot, takeupBlendDDot] = ...
            smoothStep5WithDerivatives(takeupU);
        positions(:, i) = (1 - takeupBlend) * hoverPosition ...
            + takeupBlend * takeupPosition;
        velocities(:, i) = (takeupBlendDot / takeupDuration) ...
            * (takeupPosition - hoverPosition);
        accelerations(:, i) = (takeupBlendDDot / takeupDuration^2) ...
            * (takeupPosition - hoverPosition);
    elseif strcmp(mode, 'LANDING_RELEASE')
        positions(:, i) = attach + cfg.takeoff.groundRadialOffset * radial;
        positions(3, i) = cfg.takeoff.groundZ ...
            - cfg.takeoff.vehicleGroundClearance;
    else
        positions(:, i) = attach + cfg.takeoff.groundRadialOffset * radial;
        positions(3, i) = cfg.takeoff.groundZ ...
            - cfg.takeoff.independentHoverHeight;
    end
end
end

function position = vehicleTakeupPosition(attach, radial, loadRotation, index, cfg)
% 收紧位置：从挂点沿**平衡绳向的反方向**退到 linkLength - preTensionSlack。
% q_i 由无人机指向负载 ⇒ 无人机 = 挂点 - d*q_i。
% ★ 用平衡绳向（`cfg.link.takeupLinkUnitsBody`，由控制器给出）而不是统一的
%   groundRadialOffset，这样交接瞬间实际绳向 == 期望绳向（eqi ≡ 0）。
%   qBody 为空（控制器求平衡解失败）时退回旧的统一外张几何并已发过 warning。
targetDistance = cfg.link.length - cfg.takeoff.preTensionSlack;
qBodyAll = cfg.link.takeupLinkUnitsBody;
if ~isempty(qBodyAll)
    qInertial = normalizeVector(loadRotation * qBodyAll(:, index));
    position = attach - targetDistance * qInertial;
    return
end
horizontalDistance = min(cfg.takeoff.groundRadialOffset, 0.85 * targetDistance);
verticalDistance = sqrt(max(targetDistance^2 - horizontalDistance^2, 0));
position = attach + horizontalDistance * radial;
position(3) = attach(3) - verticalDistance;
end

function qBody = equilibriumLinkUnitsBody(state, cfg)
%EQUILIBRIUMLINKUNITSBODY 求悬停平衡时控制器**期望**的绳向（负载体系 3 x n）。
%
% 做法：构造一个"悬停平衡"的 desired，让控制器算一次分配，取它的期望绳向：
%   * ex = 0、ev = 0（desired 取实测位置/零速度）⇒ Fd = -m0*g*e3；
%   * eR0 = 0、eOmega0 = 0（desired 姿态取实测）⇒ Md = 0。
%   此时分配结果只由挂点几何与零空间"外张偏置"决定，**与增益无关**
%   （kx/kv/ki/kR/kOmega 改多少都不影响这个平衡绳向）。
%
% ★ 为什么不在这里重抄一遍分配公式：分配（含外张偏置、零空间投影）是控制器的
%   内部约定，重抄一份必然随控制器改动而漂移 —— 本项目已经因为"同一个量被两处
%   定义"踩过坑（defaultLinkUnits 与 vehicleTakeupPosition 给出 7.2° 与 16.0° 两种
%   "收紧绳向"）。直接把控制器当作平衡绳向的唯一权威最省事也最不容易错。
% ★ 返回**负载体系**下的方向：分配本身就是在体系里做的（P 用 rho_i，rhs6 用
%   R0'*Fd），所以这样得到的收紧几何与"负载当时转了多少"无关，天然自洽。
n = cfg.vehicle.count;
qBody = [];
desired = struct();
desired.position = state.loadPosition;
desired.velocity = zeros(3, 1);
desired.acceleration = zeros(3, 1);
desired.rotation = state.loadRotation;
desired.bodyRate = zeros(3, 1);
desired.bodyRateDot = zeros(3, 1);
% ★★ 悬停平衡是**静态**构型 ⇒ 三阶导数与负载角加速度的导数都是 0。
%   但控制器的 analytic 前馈链（analyticForceCommandDerivative）**要求这两个字段存在**，
%   缺了就直接 error ⇒ 探针失败 ⇒ 收紧段退回 groundRadialOffset 几何 ⇒
%   交接瞬间绳向阶跃回到 30°+（本函数存在的意义就是为了消掉它）。
%   ✗ 曾经这里只填了 position/velocity/acceleration/rotation/bodyRate/bodyRateDot，
%     ⇒ 一运行就报"求悬停平衡绳向失败"。
desired.jerk = zeros(3, 1);
desired.bodyRateDDot = zeros(3, 1);
% 一次性 memory 副本：**与仿真主循环共用同一份初始化**。
% ✗ 原来这里只手填了 4 个字段，而控制器实际用到 16 个 —— 每给控制器加一个
%   memory 字段，探针就会先漏掉它（症状同上：静默退回 groundRadialOffset）。
probe = emptyControllerMemory(n, cfg);
try
    command = crazyflie_slung_controller(state, desired, probe, cfg);
    qInertial = command.desiredLinkUnits;
    if any(~isfinite(qInertial(:)))
        error('crazyflie_slung_simulation:NonFiniteEquilibrium', ...
            '悬停平衡绳向含非有限值。');
    end
    qBody = state.loadRotation.' * qInertial;
catch err
    qBody = [];
    warning('crazyflie_slung_simulation:TakeupGeometryFallback', ...
        ['求悬停平衡绳向失败（%s）⇒ 收紧段退回 groundRadialOffset 几何，' ...
         '交接瞬间会有较大的绳向阶跃（实测可达 36°）。'], err.message);
end
end

% ======================================================================
function memory = emptyControllerMemory(n, cfg)
%EMPTYCONTROLLERMEMORY 控制器跨步状态的一份**干净副本**。
%
% ★★ 为什么要抽成函数：仿真主循环和 equilibriumLinkUnitsBody 的"悬停平衡探针"
%    都需要一份全新的 memory。**探针原来只手填了 4 个字段，而控制器实际用到 16 个**，
%    于是每给控制器新增一个 memory 字段、探针就会先漏掉它：调用直接报错 ⇒
%    catch 到 ⇒ 收紧段静默退回 groundRadialOffset 几何 ⇒ 交接瞬间绳向阶跃回到 30°+
%    （而这个探针存在的唯一目的就是消掉那个阶跃）。
%    ⇒ 两边**共用同一份初始化**，以后加字段只会加一次，不会再漂移。
%
% ★ 字段清单必须与 crazyflie_slung_controller.m 里所有 `memory.*` 的用法一一对应。
%   加新字段时两边一起改；`_verify_python/_lint_matlab.py` 的"字段检查"可以发现
%   控制器读了但这里没初始化的名字。
memory = struct();
% --- 位置环与绳向环 ---
memory.positionIntegral = zeros(3, 1);
memory.linkIntegrals = zeros(3, n);
% 用初始绳向而不是 []：控制器在 mu_id 退化（norm < eps）时会回退到它，
% 给 [] 会在那种情况下直接索引越界。
memory.previousLinkUnits = cfg.initial.linkUnits;
% --- 独立起飞/降落控制器 ---
memory.independentPositionIntegral = zeros(3, n);
% --- 速率环积分器（PI 的 I 部分，见 cfg.rateLoop）---
memory.rateIntegral = zeros(3, n);
% --- 姿态外环：前馈（Omega_ic 的来源）与可选积分器，见 cfg.attitudeController ---
memory.previousVehicleRc = zeros(3, 3, n);
memory.previousVehicleU = zeros(3, n);
memory.feedforwardRate = zeros(3, n);
memory.attitudeIntegral = zeros(3, n);
memory.feedforwardStarted = false;
% --- 指令滤波与高增益观测器（analytic omegaC 的辅助状态）---
memory.filteredVehicleU = zeros(3, n);
memory.filteredVehicleUDot = zeros(3, n);
memory.commandFilterStarted = false(1, n);
memory.analyticForceCommandDot = zeros(3, n);
memory.analyticForceCommandDotValid = false;
memory.highGainVehicleU = zeros(3, n);
memory.highGainVehicleUDot = zeros(3, n);
memory.highGainObserverStarted = false(1, n);
end

function [distance, slack, qAll, qdAll] = ropeGeometryFromVehicles(state, cfg)
n = cfg.vehicle.count;
rhoAll = cfg.payload.attachPoints;
distance = zeros(1, n);
slack = zeros(1, n);
qAll = zeros(3, n);
qdAll = zeros(3, n);
for i = 1:n
    attach = state.loadPosition + state.loadRotation * rhoAll(:, i);
    attachVelocity = state.loadVelocity + state.loadRotation * ...
        cross(state.loadBodyRate, rhoAll(:, i));
    ropeVector = attach - state.vehiclePosition(:, i);
    distance(i) = norm(ropeVector);
    slack(i) = cfg.link.length - distance(i);
    if distance(i) < 1e-9
        qAll(:, i) = [0; 0; 1];
        qdAll(:, i) = zeros(3, 1);
    else
        qi = ropeVector / distance(i);
        relativeVelocity = attachVelocity - state.vehicleVelocity(:, i);
        qdi = (eye(3) - qi * qi.') * relativeVelocity / distance(i);
        qAll(:, i) = qi;
        qdAll(:, i) = qdi;
    end
end
end

function radial = horizontalRadial(rho, index, count)
radial = [rho(1); rho(2); 0];
if norm(radial) < 1e-9
    angle = 2 * pi * (index - 1) / max(count, 1);
    radial = [cos(angle); sin(angle); 0];
else
    radial = radial / norm(radial);
end
end

function code = modeCode(mode)
switch mode
    case 'SLACK'
        code = 0;
    case 'TAKEUP'
        code = 1;
    case 'TAUT_RAMP'
        code = 2;
    case 'ACTIVE'
        code = 3;
    case 'LANDING_TAUT'
        code = 4;
    case 'LANDING_RELEASE'
        code = 5;
    otherwise
        code = -1;
end
end

function value = smoothStep5(value)
value = min(max(value, 0), 1);
value = value.^3 .* (10 - 15 * value + 6 * value.^2);
end

function [value, firstDerivative, secondDerivative, thirdDerivative] = ...
    smoothStep5WithDerivatives(value)
% 五次 smoothstep 及其对无量纲参数的前三阶导数。
value = min(max(value, 0), 1);
firstDerivative = 30 * value.^2 .* (1 - value).^2;
secondDerivative = 60 * value .* (1 - value) .* (1 - 2 * value);
thirdDerivative = 60 - 360 * value + 360 * value.^2;
value = value.^3 .* (10 - 15 * value + 6 * value.^2);
end

function value = smoothStep5ThirdDerivative(value)
value = min(max(value, 0), 1);
value = 60 - 360 * value + 360 * value.^2;
end

function halfHeight = payloadGroundHalfHeight(R, sizeXYZ)
% 当前姿态下长方体沿惯性 z 轴的竖直半包络高度。
sizeXYZ = sizeXYZ(:);
halfHeight = 0.5 * sum(abs(R(3, :)).' .* sizeXYZ);
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
function value = wrapAngle(value)
% 将角度误差包到 [-pi, pi]
value = atan2(sin(value), cos(value));
end
