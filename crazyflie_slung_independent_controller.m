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
    omegaCmd = -cfg.attitudeController.kR(:) .* eR;
    omegaCmd = clampVector(omegaCmd, ...
        -cfg.takeoff.independentMaxBodyRate(:), ...
        cfg.takeoff.independentMaxBodyRate(:));
    thrust = clamp(-dot(A, R * e3), 0, cfg.vehicle.maxTotalThrust);

    forceAll(:, i) = A;
    rotationCommandAll(:, :, i) = Rc;
    bodyRateCommandAll(:, i) = omegaCmd;
    thrustAll(i) = thrust;
    attitudeErrorAll(:, i) = eR;
    positionErrorAll(:, i) = ep;
    velocityErrorAll(:, i) = ev;
end

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
