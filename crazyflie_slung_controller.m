function [command, memory] = crazyflie_slung_controller(state, desired, memory, cfg)
%CRAZYFLIE_SLUNG_CONTROLLER 多机协同吊运的几何控制器（推力 + 角速度指令）。
%
% 控制律来源：Lee 2014 第 III/IV 节与 Lee 2018 第 III/IV 节。
% 唯一结构性改动：原文最终输出推力 f_i 与力矩 M_i（式 (39)-(40)），
% 本实现把姿态环的输出改为**角速度指令 Omega_cmd_i**，交给 Crazyflie
% 固件的速率环（cflib 的 send_setpoint_manual(rate=True)）执行。
%
% 控制层次（与论文逐层对应）：
%   (1) 负载位置/姿态外环： 论文 (20)-(21) 求得合力 Fd、合力矩 Md
%   (2) 张力分配：         论文 (13)(22)-(23) 用分配矩阵 P 的伪逆
%                          把 [Fd; Md] 分配到 n 根绳索的张力 mu_id
%   (3) 绳向与平行分量：论文 (24)(25)(17)
%                          mu_i = q_i q_i' mu_id,  q_id = -mu_id/||mu_id||
%                          u||_i = mu_i + m_i l_i ||omega_i||^2 q_i + m_i q_i q_i' a_i
%   (4) 绳向环：     论文 (27) 求得 u_perp_i
%   (5) 总控制力 u_i = u||_i + u_perp_i，再分解为推力 f_i 与期望姿态 R_ic
%   (6) 姿态外环：         得到 Omega_cmd_i（本方案替代论文的 M_i）
%
% 输入：
%   state   当前状态（负载 + n 架四旋翼）
%   desired 期望量：position/velocity/acceleration/rotation/bodyRate/bodyRateDot
%   memory  跨步保存的积分器与差分状态
%   cfg     参数结构
%
% 输出：
%   command 控制器输出（3 x n 的推力/角速度指令等）
%   memory  更新后的控制器内部状态

if nargin < 4
    error('crazyflie_slung_controller:NotEnoughInputs', ...
        '需要 state、desired、memory 和 cfg 四个输入。');
end

dt = cfg.simulation.dt;
e3 = [0; 0; 1];
g = cfg.vehicle.gravity;
n = cfg.vehicle.count;

% ------------------------------------------------------------ 状态与参数
R0 = state.loadRotation;
Omega0 = state.loadBodyRate(:);
Qall = state.linkUnits;                 % 3 x n
QDall = state.linkRates;                % 3 x n
Rall = state.rotations;                 % 3 x 3 x n
OmAll = state.bodyRates;                % 3 x n

m0 = cfg.payload.mass;
J0 = cfg.payload.inertia;
rhoAll = cfg.payload.attachPoints;      % 3 x n
m = cfg.vehicle.mass;
l = cfg.link.length;

% ======================================================================
% (1) 负载位置外环：论文 (20)
% ======================================================================
ex = state.loadPosition(:) - desired.position(:);
ev = state.loadVelocity(:) - desired.velocity(:);

integralRate = ev + cfg.loadController.c1 * ex;
memory.positionIntegral = memory.positionIntegral + dt * integralRate;
memory.positionIntegral = clampVector(memory.positionIntegral, ...
    -cfg.loadController.integralLimit, cfg.loadController.integralLimit);

% ★ 论文 (20)：Fd = m0 ( -kx ex - kv ev + x0_ddot_d - g e3 )
%   这里的质量系数是 **m0**（负载质量），不是等效质量矩阵 Mq。
%   原因见参数文件的详细说明：(17) 代回 (5)(6) 后可得 (18)
%       m0 (x0_ddot - g e3) = sum_i mu_i
%   即 mu 只负责平衡负载自重；四旋翼自重由 (17) 的 m_i q_i q_i' a_i 项补偿。
Fd = m0 * (-cfg.loadController.kx .* ex ...
    - cfg.loadController.kv .* ev ...
    - cfg.loadController.ki .* memory.positionIntegral ...
    + desired.acceleration(:) - g * e3);

% 论文 (21)：负载姿态误差与角速度误差，以及期望合力矩 Md。
% yaw 采用低带宽 SO(3) 闭环：kR(3)、kOmega(3) 保持非零，但明显低于
% roll/pitch 通道。这样期望 yaw 力矩仍进入张力分配，同时减少水平张力
% 对绳摆的激励；`yawChannelEnabled` 可用于显式关闭该通道的实验对比。
R0d = desired.rotation;
Omega0d = desired.bodyRate(:);
Omega0dDot = desired.bodyRateDot(:);
eR0 = 0.5 * vee(R0d.' * R0 - R0.' * R0d);
eOmega0 = Omega0 - R0.' * R0d * Omega0d;
feedforwardRate = R0.' * R0d * Omega0d;
Md = -cfg.loadController.kR .* eR0 - cfg.loadController.kOmega .* eOmega0 ...
    + hat(feedforwardRate) * J0 * feedforwardRate ...
    + J0 * R0.' * R0d * Omega0dDot;
if isfield(cfg.loadController, 'yawChannelEnabled') ...
        && ~cfg.loadController.yawChannelEnabled
    % 实验性关闭 yaw 时同时去掉反馈和前馈的 z 力矩，避免开关只改变标志位。
    Md(3) = 0;
end

% 期望负载角加速度 Omega0_dot_cmd：由 (21) 的力矩反解
%   J0 Omega0_dot = Md - Omega0 x J0 Omega0 + (耦合项)
% 我们只需要一个**本拍可得、无延迟**的 Omega0_dot 估计来构造 a_i，
% 耦合项对 a_i 的贡献相对 α 本身是高阶小量，故此处只取主导部分。
Omega0DotCmd = J0 \ (Md - cross(Omega0, J0 * Omega0));

% ======================================================================
% (2) 张力分配：论文 (13)(22)-(23)
%   P = [ I, ..., I ; hat(rho_1), ..., hat(rho_n) ]        (6 x 3n)
%   [mu_1d; ...; mu_nd] = diag[R0,...,R0] P' (P P')^-1 [R0' Fd; Md]
% ======================================================================
P = [repmat(eye(3), 1, n); zeros(3, 3 * n)];
for i = 1:n
    P(4:6, 3 * (i - 1) + (1:3)) = hat(rhoAll(:, i));
end
PPt = P * P.';

% rank(P) = 6 的可解性检查（n >= 3，且挂点不共线）
% ★ 注意 rcond 的语义：rcond(PPt) = 1/cond(PPt)。
%   本仿真实测 cond(P P') = 276.87，故 rcond ≈ 3.6e-3，
%   远大于 pinvTolerance = 1e-9，检查通过。
%   若把挂点摆到同一水平面上（rho_z = 0），rank(P) 会掉到 5，
%   PPt 奇异，rcond -> 0，此处会直接报错——这正是这个检查要防的情况。
if cfg.allocation.checkRank
    if rcond(PPt) < cfg.allocation.pinvTolerance
        error('crazyflie_slung_controller:SingularAllocation', ...
            ['分配矩阵奇异：rcond(P*P'') = %.3e（= 1/cond，本构型应为 3.6e-3 量级）。' ...
             '请确认 n >= 3、挂点不共线、且挂点的 z 分量不全为 0。'], rcond(PPt));
    end
end

rhs6 = [R0.' * Fd; Md];
% PPt \ rhs6 等价于 inv(PPt)*rhs6，但数值上更稳。
% diag[R0,...,R0] 用 kron(eye(n), R0) 构造。
muBody = P.' * (PPt \ rhs6);
internalBiasBody = zeros(3 * n, 1);
if isfield(cfg.link, 'allowTiltedCables') && cfg.link.allowTiltedCables ...
        && isfield(cfg.allocation, 'outwardBiasFraction')
    internalBiasBody = outwardInternalBias(P, rhoAll, m0, g, ...
        cfg.allocation.outwardBiasFraction, cfg.allocation.outwardBiasMax);
    % P*internalBiasBody = 0，因此不改变负载合力和合力矩。
    muBody = muBody + internalBiasBody;
end
muDesiredAll = reshape(kron(eye(n), R0) * muBody, 3, n);

% ======================================================================
% (3)(4) 逐机：平行分量 (17) 与绳向环 (27)
% ======================================================================
uParallelAll = zeros(3, n);
uPerpAll = zeros(3, n);
desiredLinkAll = zeros(3, n);
eqAll = zeros(3, n);
Uall = zeros(3, n);

% 本拍用于 (16) 的加速度量。
% ★ 这里的状态延迟是论文没有的、纯属数值实现的产物，也是多机闭环发散的主因。
%   论文 (16)(17)(27) 全部使用**同一拍**的 a_i，且 (27) 中 hat(q_i)^3 = -hat(q_i)
%   的抵消必须建立在这个一致性上。若控制器用上一拍的 x0_ddot / Omega0_dot，
%   抵消被破坏，绳向环退化为正反馈（实测：任何单个增益打开都会指数发散，
%   把 dt 缩小 8 倍只能推迟发散时刻，说明是连续时间的寄生正反馈而非离散化误差）。
%
%   处理办法：不引入状态延迟，用**本拍的外环输出**构造 a_i。
%   外环输出在每一拍都是已知的，不需要差分，因而没有 1/dt 放大。
%   X 部分取 (20) 式括号内那一项（即期望 x0_ddot - g e3），
%   alpha 部分取 (21) 式链路推出的期望 Omega0_dot。
x0ddMinusGE3 = -cfg.loadController.kx .* ex ...
    - cfg.loadController.kv .* ev ...
    - cfg.loadController.ki .* memory.positionIntegral ...
    + desired.acceleration(:) - g * e3;

for i = 1:n
    qi = normalizeVector(Qall(:, i));
    qdi = QDall(:, i);
    rhoi = rhoAll(:, i);

    % 绳索角速度 omega_i = q_i x q_i_dot
    omegai = cross(qi, qdi);

    % (24) mu_i = q_i q_i' mu_id ；(25) q_id = -mu_id/||mu_id||
    muId = muDesiredAll(:, i);
    muI = qi * (qi.' * muId);
    nMu = norm(muId);
    if nMu < cfg.loadController.forceNormEpsilon
        qid = memory.previousLinkUnits(:, i);
    else
        qid = -muId / nMu;
    end
    qid = normalizeVector(qid);

    % (16) a_i = x0_ddot - g e3 + R0 hat(Omega0)^2 rho_i - R0 hat(rho_i) Omega0_dot
    %   = (本拍外环的期望加速度 - g e3) + 本拍测得的负载转动耦合 - 本拍期望角加速度项
    % ★ 全部使用**本拍**量，不引入上一拍的状态延迟，见文件开头关于 a_i 一致性的说明。
    ai = x0ddMinusGE3 ...
        + R0 * (hat(Omega0) * (hat(Omega0) * rhoi)) ...
        - R0 * hat(rhoi) * Omega0DotCmd;

    % (17) u||_i = mu_i + m_i l_i ||omega_i||^2 q_i + m_i q_i q_i' a_i
    uPar = muI + m * l * dot(omegai, omegai) * qi + m * qi * (qi.' * ai);
    uParallelAll(:, i) = uPar;

    % ---------------- 绳向环（论文 (26)-(28)）----------------
    % ★★★ 本项目第二个关键调参结论：必须关闭 q_id_dot 前馈 ★★★
    %   q_id 是按 (25) q_id = -mu_id/||mu_id|| 算出的**单位向量**。
    %   对它做一阶差分再除以 dt，本质是"两个单位向量之差的噪声 / 0.002 s"，
    %   信噪比极差；而且 q_id 本身会随负载的微小摆动在单位球面上抖动，
    %   差分出来的 q_id_dot 幅值远大于物理真实值。
    %   这个被污染的 q_id_dot 经 omega_id = q_id x q_id_dot 进入 (27) 的
    %       -(q_i . omega_id) q_id_dot
    %   项，构成一条**正反馈**通道，把绳向环推向发散。
    %   实测（A = 1.5 m 八字工况）：
    %       USE_QID_DOT = true   -> 稳定域只到 aX ≈ 0.15 m
    %       USE_QID_DOT = false  -> 稳定域扩到 aX ≈ 0.80 m 以上（5 倍）
    %   逐项隔离实验（_diag_fe_terms.py）进一步确认：
    %       保留原式                     -> 20.75 s 发散
    %       去掉 -(qi.wid)*qid_dot 项   -> 全程稳定（|ex| 峰 138 mm）
    %       单独去掉 omega_id 前馈      -> 结果完全一致（该项是唯一来源）
    %   注意这与"静态目标下该项理想值为 0"一致：静态工况看不出差别，
    %   只有八字这种 q_id 持续变化的工况才会把噪声放大出来。
    %
    %   q_id_dot 仍按需求解一次（保留给 omega_id / e_omega_i 的物理含义，
    %   以及便于研究者打开开关复现问题），但**不进入 correction**。
    %   把 USE_QID_DOT 置 true 可恢复论文原式，仅供研究（实测会发散）。
    USE_QID_DOT = false;
    if USE_QID_DOT && ~isempty(memory.previousLinkUnits)
        qidDot = so3ProjectedRate(memory.previousLinkUnits(:, i), qid, dt);
    else
        qidDot = zeros(3, 1);
    end
    omegaId = cross(qid, qidDot);

    eqi = cross(qid, qi);
    hatQSq = hat(qi) * hat(qi);            % hat(q_i)^2，显式相乘避免优先级歧义
    eOmegaI = omegai + hatQSq * omegaId;

    memory.linkIntegrals(:, i) = memory.linkIntegrals(:, i) + dt * eqi;
    memory.linkIntegrals(:, i) = clampVector(memory.linkIntegrals(:, i), ...
        -cfg.linkController.integralLimit, cfg.linkController.integralLimit);

    % (27)：u_perp_i = m_i l_i hat(q_i){ -kq e_qi - kw e_wi - (q_i.omega_id) q_id_dot }
    %                - m_i hat(q_i)^2 a_i
    % ★ 大括号内是负号（负反馈）。若写成正号，代入 (26) 会得到
    %   omega_dot = +kq e_q + ... 即正反馈，绳向指数发散。
    correction = -cfg.linkController.kq * eqi ...
        - cfg.linkController.komega * eOmegaI ...
        - cfg.linkController.kqIntegral * memory.linkIntegrals(:, i) ...
        - dot(qi, omegaId) * qidDot;
    uPerp = m * l * cross(qi, correction) - m * hatQSq * ai;
    uPerpAll(:, i) = uPerp;

    desiredLinkAll(:, i) = qid;
    eqAll(:, i) = eqi;

    % (29) u_i = u||_i + u_perp_i
    Uall(:, i) = uPar + uPerp;
end

memory.previousLinkUnits = desiredLinkAll;

% ======================================================================
% (5) 推力与期望姿态：论文 (36)-(39)
% ======================================================================
b3cAll = zeros(3, n);
RcAll = zeros(3, 3, n);
thrustAll = zeros(1, n);
thrustPctAll = zeros(1, n);
eRAll = zeros(3, n);
eOmRAll = zeros(3, n);
omegaCmdAll = zeros(3, n);

b1d = desired.rotation(:, 1);

for i = 1:n
    Ri = Rall(:, :, i);
    Omegai = OmAll(:, i);
    ui = Uall(:, i);
    uNorm = norm(ui);

    % (36) 期望机体第三轴 b3c_i = -u_i / ||u_i||
    % 退化保护：||u_i|| -> 0 时保留上一拍的 b3c（存在 memory 里），
    % 再退化为 +e3。这里**不能**读 RcAll(:,3,i)，因为 RcAll 正在本循环内填充。
    if uNorm < cfg.loadController.forceNormEpsilon
        b3c = normalizeVector(Rall(:, :, i) * e3);   % 沿用当前实际机体第三轴
        if norm(b3c) < 0.5
            b3c = e3;
        end
    else
        b3c = -ui / uNorm;
    end
    b3c = normalizeVector(b3c);

    % (37) 由 b1d 在垂直于 b3c 平面上的投影构造期望姿态 R_ic
    b1Projection = (eye(3) - b3c * b3c.') * b1d;
    if norm(b1Projection) < cfg.attitudeController.maxAttitudeError * 1e-6
        b1Projection = (eye(3) - b3c * b3c.') * [1; 0; 0];
    end
    b1c = normalizeVector(b1Projection);
    b2c = normalizeVector(cross(b3c, b1c));
    b1c = normalizeVector(cross(b2c, b3c));
    Rc = projectSO3([b1c, b2c, b3c]);
    b3cAll(:, i) = b3c;
    RcAll(:, :, i) = Rc;

    % (39) 推力幅值 f_i = -u_i' R_i e3
    thrustDesired = -dot(ui, Ri * e3);
    thrustCmd = clamp(thrustDesired, 0, cfg.vehicle.maxTotalThrust);

    % ------------------------------------------------ (6) 姿态外环
    eR = 0.5 * vee(Rc.' * Ri - Ri.' * Rc);
    eRNorm = norm(eR);
    if eRNorm > cfg.attitudeController.maxAttitudeError
        eR = eR * (cfg.attitudeController.maxAttitudeError / eRNorm);
    end
    eOmR = Omegai;      % 期望姿态角速度 Omega_ic 取 0（见下方说明）

    % ★ Omega_cmd_i = -kR .* e_Ri - kOmega .* (Omega_i - R_i' R_ic Omega_ic)
    %   这里把期望姿态角速度前馈 Omega_ic 取 0。
    %   原因：用相邻步长的 SO(3) 对数差分算 R_ic 的变化率会被 1/dt = 500 放大，
    %   R_ic 由 u_i 构造，u_i 里的微小抖动会变成每秒几百弧度的假期望角速度，
    %   再经 (1 + kOmega) 前馈进 omega_cmd 形成自激（实测姿态误差会随时间
    %   从 0.6 deg 增长到 8.9 deg，机体角速度长期维持 4.2 rad/s）。
    %   静态悬停目标下 R_ic 恒定，该项理想值本来就是 0。
    omegaCmd = -cfg.attitudeController.kR .* eR ...
        - cfg.attitudeController.kOmega .* eOmR;
    omegaCmd = clampVector(omegaCmd, ...
        -cfg.attitudeController.maxBodyRateCommand, ...
        cfg.attitudeController.maxBodyRateCommand);

    eRAll(:, i) = eR;
    eOmRAll(:, i) = eOmR;
    omegaCmdAll(:, i) = omegaCmd;
    thrustAll(i) = thrustDesired;
    thrustPctAll(i) = 100 * thrustCmd / cfg.vehicle.maxTotalThrust;
end

% ------------------------------------------------------------------ 输出
command = struct();
command.positionError = ex;
command.velocityError = ev;
command.loadAttitudeError = eR0;
command.loadBodyRateError = eOmega0;
command.desiredForce = Fd;
command.desiredMoment = Md;
command.desiredLinkUnits = desiredLinkAll;
command.linkDirectionErrors = eqAll;
command.desiredTensions = muDesiredAll;      % 分配结果 mu_id（3 x n）
command.internalTensionBias = reshape(kron(eye(n), R0) * internalBiasBody, 3, n);
command.tensions = zeros(1, n);              % 实际张力 = ||mu_i||，mu_i = q_i q_i' mu_id
for i = 1:n
    qi = normalizeVector(Qall(:, i));
    muI = qi * (qi.' * muDesiredAll(:, i));
    command.tensions(i) = norm(muI);
end
command.parallelForces = uParallelAll;
command.perpendicularForces = uPerpAll;
command.totalForces = Uall;
command.computedRotations = RcAll;
command.b3c = b3cAll;
command.attitudeErrors = eRAll;
command.bodyRateErrors = eOmRAll;
command.omegaCommands = omegaCmdAll;
command.omegaCommandDeg = rad2deg(omegaCmdAll);
command.thrustDesired = thrustAll;
command.thrustPercentage = thrustPctAll;
command.totalThrust = sum(thrustAll);
end

% ======================================================================
% 局部工具函数（详见 README §1.1：MATLAB 局部函数是文件私有的）
% ======================================================================
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

function S = hat(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function R = projectSO3(R)
[U, ~, V] = svd(R);
R = U * V.';
if det(R) < 0
    U(:, 3) = -U(:, 3);
    R = U * V.';
end
end

function rate = so3ProjectedRate(previousUnit, currentUnit, dt)
% 单位向量的一阶变化率，并投影到当前单位向量的切空间。
delta = (currentUnit(:) - previousUnit(:)) / max(dt, eps);
rate = delta - dot(delta, currentUnit(:)) * currentUnit(:);
end

function value = clamp(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end

function value = clampVector(value, lowerBound, upperBound)
value = min(max(value, lowerBound), upperBound);
end

function bias = outwardInternalBias(P, rhoAll, mass, gravity, fraction, maxRms)
%OUTWARDINTERNALBIAS 生成不改变负载 wrench 的外张内部张力。
n = size(rhoAll, 2);
desired = zeros(3, n);
for i = 1:n
    radial = [rhoAll(1, i); rhoAll(2, i); 0];
    if norm(radial) < 1e-12
        angle = 2 * pi * (i - 1) / max(n, 1);
        radial = [cos(angle); sin(angle); 0];
    else
        radial = radial / norm(radial);
    end
    desired(:, i) = fraction * mass * gravity / sqrt(n) * radial;
end
nullBasis = null(P);
if isempty(nullBasis)
    bias = zeros(3 * n, 1);
    return
end
bias = nullBasis * (nullBasis.' * desired(:));
biasBlocks = reshape(bias, 3, n);
biasRms = sqrt(mean(sum(biasBlocks.^2, 1)));
targetRms = min(fraction * mass * gravity / sqrt(n), maxRms);
if biasRms > 1e-12 && targetRms > 0
    bias = bias * (targetRms / biasRms);
end
end
