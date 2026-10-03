function [command, memory] = crazyflie_slung_independent_controller(state, desiredPosition, desiredVelocity, desiredAcceleration, memory, cfg)
%CRAZYFLIE_SLUNG_INDEPENDENT_CONTROLLER
% 起飞/释放阶段的单机几何位置控制器。
%
% 该控制器只在绳索松弛时使用。每架无人机独立跟踪一个安全位置，采用
% v2 中的几何 PID -> 期望推力方向 -> SO(3) 姿态误差 -> 角速度指令链路。
% 绳索松弛时不把负载期望力分配给无人机，也不调用绷紧段的张力控制器。

if nargin < 6
    error('crazyflie_slung_independent_controller:NotEnoughInputs', ...
        '需要 state、desiredPosition、desiredVelocity、desiredAcceleration、memory 和 cfg。');
end

dt = cfg.simulation.dt;
e3 = [0; 0; 1];
n = cfg.vehicle.count;
m = cfg.vehicle.mass;
g = cfg.vehicle.gravity;
Kp = cfg.takeoff.independentPositionKp(:);
Kv = cfg.takeoff.independentPositionKv(:);
if isfield(cfg.takeoff, 'independentIntegralGain')
    Ki = cfg.takeoff.independentIntegralGain(:);
else
    Ki = 0.8 * ones(3, 1);
end
heading = normalizeVector(cfg.takeoff.independentHeading(:));

if ~isfield(memory, 'independentPositionIntegral') || ...
        ~isequal(size(memory.independentPositionIntegral), [3, n])
    memory.independentPositionIntegral = zeros(3, n);
end

desiredPosition = reshape(desiredPosition, 3, n);
desiredVelocity = reshape(desiredVelocity, 3, n);
desiredAcceleration = reshape(desiredAcceleration, 3, n);

forceAll = zeros(3, n);
rotationCommandAll = zeros(3, 3, n);
bodyRateCommandAll = zeros(3, n);
thrustAll = zeros(1, n);
attitudeErrorAll = zeros(3, n);
positionErrorAll = zeros(3, n);
velocityErrorAll = zeros(3, n);

for i = 1:n
    p = state.vehiclePosition(:, i);
    v = state.vehicleVelocity(:, i);
    R = projectSO3(state.rotations(:, :, i));

    ep = p - desiredPosition(:, i);
    ev = v - desiredVelocity(:, i);
    integralRate = ev + 0.5 * ep;
    % ★★ 抗积分饱和（2026-09-29 新增）：位置误差大时**不积分**。
    %   原来无条件积分，而 SLACK/TAKEUP 是"从 0.3 m 外飞向目标"的大机动：
    %   一路灌到限幅附近（Ki=0.8、限幅 0.20 ⇒ 最大 0.16 m/s² 的恒定力偏置），
    %   等停下来时这个偏置还在 ⇒ 无人机被稳稳压在目标**下方 18 mm**。
    %   实测后果：收紧结束时绳长停在 0.6200 m = l − 0.0300，**正好落在
    %   epsilonOn = 0.030 的判据边界上** ⇒ "绳是否绷紧"变成临界判断，
    %   确认时间被拖长（TAKEUP 实际持续 3.28 s 而不是 takeupDuration 的 2.5 s）。
    %   加上 gate 后：大机动段积分器不动（不影响跟踪），
    %   误差进入 gate 后才积分 ⇒ 稳态余差被收干净，收紧长度能到设计值 l−0.012。
    if norm(ep) < cfg.takeoff.independentIntegralGate
        memory.independentPositionIntegral(:, i) = ...
            memory.independentPositionIntegral(:, i) + dt * integralRate;
        memory.independentPositionIntegral(:, i) = clampVector(...
            memory.independentPositionIntegral(:, i), ...
            -cfg.takeoff.independentIntegralLimit(:), ...
            cfg.takeoff.independentIntegralLimit(:));
    end

    % Kp/Kv/Ki are acceleration gains.  Convert the complete translational
    % command to force here; omitting m makes the feedback roughly 30 times
    % too large for a Crazyflie and causes the takeoff trajectory to diverge.
    feedbackAcceleration = -Kp .* ep - Kv .* ev ...
        - Ki .* memory.independentPositionIntegral(:, i);
    if isfield(cfg.takeoff, 'independentMaxFeedbackAcceleration')
        feedbackAcceleration = clampVector(feedbackAcceleration, ...
            -cfg.takeoff.independentMaxFeedbackAcceleration(:), ...
            cfg.takeoff.independentMaxFeedbackAcceleration(:));
    end
    A = m * (feedbackAcceleration + desiredAcceleration(:, i) - g * e3);
    if isfield(cfg.attitudeController, 'omegaCMethod') ...
            && strcmpi(cfg.attitudeController.omegaCMethod, 'command_filter')
        [A, memory] = filterAttitudeForceCommand(A, i, memory, cfg);
    elseif isfield(cfg.attitudeController, 'omegaCMethod') ...
            && strcmpi(cfg.attitudeController.omegaCMethod, 'high_gain_observer')
        [~, memory] = updateHighGainForceDerivative(A, i, memory, cfg);
    end
    if norm(A) < 1e-8
        b3c = e3;
    else
        b3c = -A / norm(A);
    end

    b1Projection = (eye(3) - b3c * b3c.') * heading;
    if norm(b1Projection) < 1e-6
        b1Projection = (eye(3) - b3c * b3c.') * [1; 0; 0];
    end
    b1c = normalizeVector(b1Projection);
    b2c = normalizeVector(cross(b3c, b1c));
    b1c = normalizeVector(cross(b2c, b3c));
    Rc = projectSO3([b1c, b2c, b3c]);

    eR = 0.5 * vee(Rc.' * R - R.' * Rc);
    % 期望姿态角速度前馈（与 v2 同构，与主控制器 crazyflie_slung_controller.m 一致）：
    %   Omega_cmd = R'R_c Omega_c - kR .* e_R          （useRateDamping=false，同 v2）
    % Omega_c 由 omegaCMethod 选择：'analytic'、'filtered_log_difference' 或
    % 'command_filter'；关掉用 'none'。heading 是常量 ⇒ b1dDot = 0。
    omegaCi = feedforwardBodyRate(Rc, A, heading, zeros(3, 1), i, memory, cfg);
    feedforward = R.' * Rc * omegaCi;
    eOmR = state.bodyRates(:, i) - feedforward;
    omegaCmd = feedforward - cfg.attitudeController.kR(:) .* eR;
    if cfg.attitudeController.useRateDamping
        omegaCmd = omegaCmd - cfg.attitudeController.kOmega(:) .* eOmR;
    end
    omegaCmd = clampVector(omegaCmd, ...
        -cfg.takeoff.independentMaxBodyRate(:), ...
        cfg.takeoff.independentMaxBodyRate(:));
    memory.previousVehicleRc(:, :, i) = Rc;
    memory.previousVehicleU(:, i) = A;
    memory.feedforwardRate(:, i) = omegaCi;
    thrust = clamp(-dot(A, R * e3), 0, cfg.vehicle.maxTotalThrust);

    forceAll(:, i) = A;
    rotationCommandAll(:, :, i) = Rc;
    bodyRateCommandAll(:, i) = omegaCmd;
    thrustAll(i) = thrust;
    attitudeErrorAll(:, i) = eR;
    positionErrorAll(:, i) = ep;
    velocityErrorAll(:, i) = ev;
end
memory.feedforwardStarted = true;       % 下一拍才有可用的历史 R_ic

command = struct();
command.positionError = mean(positionErrorAll, 2);
command.velocityError = mean(velocityErrorAll, 2);
command.loadAttitudeError = zeros(3, 1);
command.loadBodyRateError = zeros(3, 1);
command.desiredForce = mean(forceAll, 2);
command.desiredMoment = zeros(3, 1);
command.desiredLinkUnits = state.linkUnits;
command.linkDirectionErrors = zeros(3, n);
command.desiredTensions = zeros(3, n);
command.internalTensionBias = zeros(3, n);
command.tensions = zeros(1, n);
command.parallelForces = zeros(3, n);
command.perpendicularForces = zeros(3, n);
command.totalForces = forceAll;
command.computedRotations = rotationCommandAll;
command.b3c = rotationCommandAll(:, 3, :);
command.attitudeErrors = attitudeErrorAll;
command.bodyRateErrors = zeros(3, n);
command.omegaCommands = bodyRateCommandAll;
command.omegaCommandDeg = rad2deg(bodyRateCommandAll);
command.thrustDesired = thrustAll;
command.thrustPercentage = 100 * thrustAll / cfg.vehicle.maxTotalThrust;
command.totalThrust = sum(thrustAll);
command.independentMode = true;
end

function q = normalizeVector(q)
q = q(:);
n = norm(q);
if n < eps
    q = [0; 0; 1];
else
    q = q / n;
end
end

% ======================================================================
function omegaC = feedforwardBodyRate(Rc, ui, b1d, b1dDot, index, memory, cfg)
%FEEDFORWARDBODYRATE 期望姿态 R_ic 的角速度前馈 Omega_ic（负载体系）。
% 与 crazyflie_slung_controller.m 里的同名函数**完全一致**（MATLAB 的局部函数是
% 文件私有的，按本工程约定各自复制一份）。三种算法见那边的说明。
omegaC = zeros(3, 1);
method = 'none';
if isfield(cfg.attitudeController, 'omegaCMethod')
    method = cfg.attitudeController.omegaCMethod;
end
if strcmp(method, 'none') ...
        || ~isfield(memory, 'feedforwardStarted') || ~memory.feedforwardStarted
    return
end

switch method
    case 'analytic'
        uiDot = (ui - memory.previousVehicleU(:, index)) / max(cfg.simulation.dt, eps);
        raw = analyticOmegaC(Rc, ui, uiDot, b1d, b1dDot, cfg);
    case 'command_filter'
        raw = analyticOmegaC(Rc, ui, ...
            memory.filteredVehicleUDot(:, index), b1d, b1dDot, cfg);
    case 'high_gain_observer'
        raw = analyticOmegaC(Rc, ui, ...
            memory.highGainVehicleUDot(:, index), b1d, b1dDot, cfg);
    case 'filtered_log_difference'
        previousRc = memory.previousVehicleRc(:, :, index);
        raw = so3Log(previousRc.' * Rc) / max(cfg.simulation.dt, eps);
    case 'none'
        return
    otherwise
        return
end

raw = clampVector(raw, -cfg.attitudeController.feedforwardMaxRate, ...
    cfg.attitudeController.feedforwardMaxRate);
tau = cfg.attitudeController.feedforwardFilterTime;
if tau > 0
    alpha = cfg.simulation.dt / (tau + cfg.simulation.dt);
else
    alpha = 1;
end
omegaC = memory.feedforwardRate(:, index) ...
    + alpha * (raw - memory.feedforwardRate(:, index));
end

function [filteredU, memory] = filterAttitudeForceCommand(commandU, index, memory, cfg)
% 二阶命令滤波器：ü_f + 2*zeta*wn*u̇_f + wn^2*u_f = wn^2*u。
dt = cfg.simulation.dt;
wn = cfg.attitudeController.commandFilterNaturalFrequency;
zeta = cfg.attitudeController.commandFilterDampingRatio;
if ~memory.commandFilterStarted(index)
    memory.filteredVehicleU(:, index) = commandU;
    memory.filteredVehicleUDot(:, index) = zeros(3, 1);
    memory.commandFilterStarted(index) = true;
else
    uFiltered = memory.filteredVehicleU(:, index);
    uFilteredDot = memory.filteredVehicleUDot(:, index);
    uFilteredDDot = wn^2 * (commandU - uFiltered) - 2 * zeta * wn * uFilteredDot;
    uFilteredDot = uFilteredDot + dt * uFilteredDDot;
    uFiltered = uFiltered + dt * uFilteredDot;
    memory.filteredVehicleU(:, index) = uFiltered;
    memory.filteredVehicleUDot(:, index) = uFilteredDot;
end
filteredU = memory.filteredVehicleU(:, index);
end

function [estimatedU, memory] = updateHighGainForceDerivative(commandU, index, memory, cfg)
% 二阶高增益非线性微分器，第二状态 z1 估计 commandU 的导数。
dt = cfg.simulation.dt;
lambda1 = cfg.attitudeController.highGainDifferentiatorLambda1;
lambda2 = cfg.attitudeController.highGainDifferentiatorLambda2;
sigma = max(cfg.attitudeController.highGainDifferentiatorSmoothing, 1e-9);
limit = cfg.attitudeController.forceDerivativeLimit;
if ~memory.highGainObserverStarted(index)
    memory.highGainVehicleU(:, index) = commandU;
    memory.highGainVehicleUDot(:, index) = zeros(3, 1);
    memory.highGainObserverStarted(index) = true;
else
    z0 = memory.highGainVehicleU(:, index);
    z1 = memory.highGainVehicleUDot(:, index);
    error = commandU - z0;
    injection = tanh(error / sigma);
    z0Dot = z1 + lambda1 * sqrt(abs(error) + 1e-12) .* injection;
    z1Dot = lambda2 * injection;
    z1 = clampVector(z1 + dt * z1Dot, -limit, limit);
    z0 = z0 + dt * z0Dot;
    memory.highGainVehicleU(:, index) = z0;
    memory.highGainVehicleUDot(:, index) = z1;
end
estimatedU = memory.highGainVehicleU(:, index);
end

% ======================================================================
function OmegaC = analyticOmegaC(Rc, ui, uiDot, b1d, b1dDot, cfg)
%ANALYTICOMEGAC 由 R_c 的列向量导数解析求 Omega_c = vee(R_c'RcDot)（与 v2 同构）。
% 独立起飞控制器没有绷紧段张力分配链；analytic 模式在此仍使用
% 独立力指令的离散导数，完整解析链只用于主协同搬运控制器。
b3c = Rc(:, 3);
forceNorm = norm(ui);
if forceNorm < cfg.loadController.forceNormEpsilon
    OmegaC = zeros(3, 1);
    return
end
b3cDot = -(eye(3) - b3c * b3c.') * uiDot / forceNorm;

b1c = Rc(:, 1);
projection = (eye(3) - b3c * b3c.') * b1d;
if norm(projection) < cfg.loadController.forceNormEpsilon
    b1cDot = zeros(3, 1);
else
    projectionDot = -(b3cDot * b3c.' + b3c * b3cDot.') * b1d ...
        + (eye(3) - b3c * b3c.') * b1dDot;
    b1cDot = (eye(3) - b1c * b1c.') * projectionDot / norm(projection);
end

b2c = Rc(:, 2);
b2cDot = cross(b3cDot, b1c) + cross(b3c, b1cDot);
RcDot = [b1cDot, b2cDot, b3cDot];
OmegaHat = Rc.' * RcDot;
OmegaHat = 0.5 * (OmegaHat - OmegaHat.');
OmegaC = vee(OmegaHat);
end

% ======================================================================
function v = so3Log(R)
%SO3LOG SO(3) 对数映射：旋转矩阵 -> 旋转向量。与 v2 的实现一致。
R = projectSO3(R);
cosTheta = min(max((trace(R) - 1) / 2, -1), 1);
theta = acos(cosTheta);
if theta < 1e-7
    v = 0.5 * vee(R - R.');
elseif abs(pi - theta) < 1e-5
    [V, D] = eig((R + eye(3)) / 2);
    [~, idx] = max(real(diag(D)));
    axisVector = real(V(:, idx));
    axisVector = axisVector / max(norm(axisVector), eps);
    v = theta * axisVector;
else
    v = theta / (2 * sin(theta)) * vee(R - R.');
end
end

function v = vee(S)
v = [S(3, 2); S(1, 3); S(2, 1)];
end

function R = projectSO3(R)
[U, ~, V] = svd(R);
R = U * V.';
if det(R) < 0
    U(:, 3) = -U(:, 3);
    R = U * V.';
end
end

function value = clamp(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end

function value = clampVector(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end
