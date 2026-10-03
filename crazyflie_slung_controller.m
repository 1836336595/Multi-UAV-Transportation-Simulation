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

% 论文 (21)：负载姿态误差与角速度误差，以及期望合力矩 Md
% 负载姿态误差采用 SO(3) 几何误差。yaw 是否参与控制由
% cfg.loadController.yawChannelEnabled 显式决定。当前默认配置开启低带宽
% yaw 反馈；Md(3) 通过张力分配和绳向环产生负载偏航恢复力矩。
R0d = desired.rotation;
Omega0d = desired.bodyRate(:);
Omega0dDot = desired.bodyRateDot(:);
eR0 = 0.5 * vee(R0d.' * R0 - R0.' * R0d);
eOmega0 = Omega0 - R0.' * R0d * Omega0d;
feedforwardRate = R0.' * R0d * Omega0d;
MdRaw = -cfg.loadController.kR .* eR0 - cfg.loadController.kOmega .* eOmega0 ...
    + hat(feedforwardRate) * J0 * feedforwardRate ...
    + J0 * R0.' * R0d * Omega0dDot;

Md = MdRaw;
if ~isfield(cfg.loadController, 'yawChannelEnabled') || ...
        ~cfg.loadController.yawChannelEnabled
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
%   挂点是否处于同一水平面并不决定秩；关键是 n >= 3 且挂点不共线。
%   只有挂点共线时 rank(P) 才会下降，PPt 才可能奇异。
if cfg.allocation.checkRank
    if rcond(PPt) < cfg.allocation.pinvTolerance
        error('crazyflie_slung_controller:SingularAllocation', ...
            ['分配矩阵奇异：rcond(P*P'') = %.3e（= 1/cond，本构型应为 3.6e-3 量级）。' ...
             '请确认 n >= 3 且挂点不共线。'], rcond(PPt));
    end
end

rhs6 = [R0.' * Fd; Md];
% PPt \ rhs6 等价于 inv(PPt)*rhs6，但数值上更稳。
% diag[R0,...,R0] 用 kron(eye(n), R0) 构造。
muBody = P.' * (PPt \ rhs6);

% 最小范数解只负责满足总合力/总力矩，悬停时会让绳索尽量接近竖直。
% 对小尺寸负载，这会让无人机中心过度聚拢。利用 P 的零空间加入内部
% “外张”力：P * muInternal = 0，因此不改变负载需要的 Fd/Md，只改变
% 各根绳的方向和无人机间距。该项是可配置的，置 0 即恢复论文的最小范数解。
internalBiasBody = zeros(3 * n, 1);
if isfield(cfg.link, 'allowTiltedCables') && cfg.link.allowTiltedCables ...
        && isfield(cfg.allocation, 'outwardBiasFraction')
    internalBiasBody = outwardInternalBias(P, rhoAll, m0, g, ...
        cfg.allocation.outwardBiasFraction, cfg.allocation.outwardBiasMax);
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
aiAll = zeros(3, n);
qidAll = zeros(3, n);
omegaAll = zeros(3, n);
correctionAll = zeros(3, n);

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
    aiAll(:, i) = ai;
    qidAll(:, i) = qid;
    omegaAll(:, i) = omegai;

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
    correctionAll(:, i) = correction;
    uPerp = m * l * cross(qi, correction) - m * hatQSq * ai;
    uPerpAll(:, i) = uPerp;

    desiredLinkAll(:, i) = qid;
    eqAll(:, i) = eqi;

    % (29) u_i = u||_i + u_perp_i
    Uall(:, i) = uPar + uPerp;
end

memory.previousLinkUnits = desiredLinkAll;

% 本工程的 analytic 前馈使用完整链式导数：
% Fd/Md -> mu_id -> q_id/mu_i -> a_i -> u_parallel/u_perp -> u_cmd。
if strcmpi(cfg.attitudeController.omegaCMethod, 'analytic')
    memory.analyticForceCommandDot = analyticForceCommandDerivative(...
        state, desired, cfg, P, PPt, muBody, muDesiredAll, Fd, Md, ...
        ex, ev, integralRate, memory.positionIntegral, aiAll, qidAll, omegaAll, correctionAll, ...
        uPerpAll, memory.linkIntegrals);
    memory.analyticForceCommandDotValid = true;
else
    memory.analyticForceCommandDotValid = false;
end

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

% (37) 期望姿态 R_ic 的第一轴参考方向 b1d。
%   'reference' = 论文原式（取 R0d 第一轴）；'worldX' = 机体航向锁定世界 +x。
% ★ 机体航向跟着负载参考 yaw 转时，挂点方位也随之变化；因此动态参考下
%   需要保证绳向环和负载 yaw 环带宽匹配。四旋翼自身 yaw 与推力方向解耦，
%   这里的 headingSource 只选择机体绕推力轴的姿态参考，不会关闭负载 yaw 环。
%   实测（八字）：偏航率抖幅 std 0.0498 → 0.0269 rad/s，终端残差 94.0 → 90.0 mm。
% ★ 定高工况下 R0d ≡ cfg.target.R0（第一轴 = +x），两者**完全等价**。
if isfield(cfg.attitudeController, 'headingSource') ...
        && strcmpi(cfg.attitudeController.headingSource, 'worldX')
    b1d = [1; 0; 0];
else
    b1d = desired.rotation(:, 1);
end

% b1d 的导数 —— 解析前馈（v2 的 analyticOmegaC）需要它。
% 本工程的 desired 只提供 bodyRate（负载期望角速度 Omega0d，机体系），
% 于是 R_dot = R0d * hat(Omega0d) ⇒ b1dDot = R0d * (Omega0d × e1)。
% 静态参考下 Omega0d = 0 ⇒ b1dDot = 0。
% ★ 若日后某个参考函数额外提供 desired.rotationDot，可在此优先取它
%   （v2 就是这么做的）；当前管线不产生该字段，所以不写死在这里，
%   否则静态检查会报"desired 没有这个字段"。
if isfield(desired, 'bodyRate') && ~isempty(desired.bodyRate)
    b1dDot = desired.rotation * cross(desired.bodyRate(:), [1; 0; 0]);
else
    b1dDot = zeros(3, 1);
end

for i = 1:n
    Ri = Rall(:, :, i);
    Omegai = OmAll(:, i);
    ui = Uall(:, i);
    if isfield(cfg.attitudeController, 'omegaCMethod') ...
            && strcmpi(cfg.attitudeController.omegaCMethod, 'command_filter')
        uiRaw = ui;
        [ui, memory] = filterAttitudeForceCommand(ui, i, memory, cfg);
        forceDelta = ui - uiRaw;
        qi = normalizeVector(Qall(:, i));
        uParallelAll(:, i) = uParallelAll(:, i) + qi * (qi.' * forceDelta);
        uPerpAll(:, i) = uPerpAll(:, i) ...
            + (eye(3) - qi * qi.') * forceDelta;
    elseif isfield(cfg.attitudeController, 'omegaCMethod') ...
            && strcmpi(cfg.attitudeController.omegaCMethod, 'high_gain_observer')
        [~, memory] = updateHighGainForceDerivative(ui, i, memory, cfg);
    end
    Uall(:, i) = ui;
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

    % ---- 期望姿态角速度前馈 Omega_ic（v2 的结构：前馈 + 比例 + 可选积分）----
    % b1d / b1dDot 是 R_c 第一轴参考及其导数，解析前馈要用
    omegaCi = feedforwardBodyRate(Rc, ui, b1d, b1dDot, i, memory, cfg);
    feedforward = Ri.' * Rc * omegaCi;      % 换到本机体系
    eOmR = Omegai - feedforward;            % 论文 (30) 的 e_Omega_i

    % ---- 外层姿态积分（对应 v2 的 useIntegral，默认同样关闭）----
    if cfg.attitudeController.useIntegral
        memory.attitudeIntegral(:, i) = memory.attitudeIntegral(:, i) + dt * eR;
        memory.attitudeIntegral(:, i) = clampVector(...
            memory.attitudeIntegral(:, i), ...
            -cfg.attitudeController.integralLimit, ...
            cfg.attitudeController.integralLimit);
    else
        memory.attitudeIntegral(:, i) = zeros(3, 1);
    end

    % ★★ useRateDamping = false（默认）= **与 v2 完全一致的式子**：
    %        Omega_cmd_i = R_i'R_ic Omega_ic - kR .* e_Ri - kIR .* ∫e_Ri
    %    v2 没有速率阻尼项（阻尼交给内环 PI）。置 true 则额外加论文 (30) 的
    %        - kOmega .* (Omega_i - R_i'R_ic Omega_ic)
    omegaCmd = feedforward ...
        - cfg.attitudeController.kR .* eR ...
        - cfg.attitudeController.kIR .* memory.attitudeIntegral(:, i);
    if cfg.attitudeController.useRateDamping
        omegaCmd = omegaCmd - cfg.attitudeController.kOmega .* eOmR;
    end
    omegaCmd = clampVector(omegaCmd, ...
        -cfg.attitudeController.maxBodyRateCommand, ...
        cfg.attitudeController.maxBodyRateCommand);

    % 记录本拍期望姿态与期望合力，供下一步算前馈（放在最后，避免自比较）
    memory.previousVehicleRc(:, :, i) = Rc;
    memory.previousVehicleU(:, i) = ui;
    memory.feedforwardRate(:, i) = omegaCi;

    eRAll(:, i) = eR;
    eOmRAll(:, i) = eOmR;
    omegaCmdAll(:, i) = omegaCmd;
    thrustAll(i) = thrustDesired;
    thrustPctAll(i) = 100 * thrustCmd / cfg.vehicle.maxTotalThrust;
end
memory.feedforwardStarted = true;       % 下一拍才有可用的历史 R_ic

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

% ======================================================================
function omegaC = feedforwardBodyRate(Rc, ui, b1d, b1dDot, index, memory, cfg)
%FEEDFORWARDBODYRATE 期望姿态 R_ic 的角速度前馈 Omega_ic（**负载体系**）。
%
% 前馈算法由 omegaCMethod 选择：
%   'analytic'（默认）        —— 按完整控制链解析传播 u_i^cmd 的导数
%   'filtered_log_difference' —— 相邻期望姿态的 SO(3) 对数差分
%   'command_filter'          —— 二阶滤波 u_i^cmd，并用滤波器状态给出导数
%   'none'                    —— 关闭前馈（Omega_ic ≡ 0）
%
% 各方法估计值经过限幅 + 一阶低通；command_filter 已先平滑力指令，
% 其角速度估计也沿用同样的限幅与低通保护。
%
% ★ 为什么仍要限幅/低通：前馈的种子（u_i^cmd 的变化 / R_c 的变化）含张力分配与
%   绳向环的高频抖动，裸的 1/dt 差分会把它放大成每秒几百弧度的假前馈（1/dt = 500）。
%   本项目**实测过这个自激**：姿态误差随时间从 0.6 deg 涨到 8.9 deg、
%   机体角速度长期维持 4.2 rad/s —— 所以历史上这里把前馈整项置零。
%   现在按 v2 结构启用，但必须低通；若仍见误差缓慢增长 ⇒ 调大
%   feedforwardFilterTime，或设 omegaCMethod = 'none'。
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
        if isfield(memory, 'analyticForceCommandDotValid') ...
                && memory.analyticForceCommandDotValid
            uiDot = memory.analyticForceCommandDot(:, index);
        else
            error('crazyflie_slung_controller:AnalyticDerivativeUnavailable', ...
                ['omegaCMethod=analytic 但本拍没有有效的解析 u_i^{cmd} 导数。' ...
                 '请检查 desired.jerk、desired.bodyRateDDot 以及解析控制链输入。']);
        end
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

% 限幅 + 一阶低通（各方法共用）
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
% 二阶高增益非线性微分器：第二状态 z1 估计 commandU 的导数。
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

function forceCommandDot = analyticForceCommandDerivative(state, desired, cfg, ...
    P, PPt, muBody, muDesiredAll, desiredForce, desiredMoment, ...
    positionError, velocityError, positionIntegralRate, positionIntegral, accelerationAll, ...
    desiredLinkAll, linkAngularVelocityAll, correctionAll, ...
    perpendicularForceAll, linkIntegralAll)
% 对 u_i = u_parallel_i + u_perp_i 逐层应用乘积法则和链式法则。
if ~isfield(desired, 'jerk') || numel(desired.jerk) ~= 3 ...
        || ~isfield(desired, 'bodyRateDDot') || isempty(desired.bodyRateDDot) ...
        || numel(desired.bodyRateDDot) ~= 3
    error('crazyflie_slung_controller:MissingReferenceDerivatives', ...
        ['analytic omegaCMethod 需要 desired.jerk (3x1) 与 ' ...
         'desired.bodyRateDDot (3x1)。请由参考轨迹解析提供这两个量；' ...
         '当前参考在非锁 yaw 巡航时尚未提供 bodyRateDDot。']);
end

n = cfg.vehicle.count;
m0 = cfg.payload.mass;
m = cfg.vehicle.mass;
l = cfg.link.length;
J0 = cfg.payload.inertia;
R0 = state.loadRotation;
Omega0 = state.loadBodyRate(:);
Omega0Dot = state.loadBodyAcceleration(:);
R0d = desired.rotation;
Omega0d = desired.bodyRate(:);
Omega0dDot = desired.bodyRateDot(:);
Omega0dDDot = desired.bodyRateDDot(:);
e3 = [0; 0; 1];

positionErrorDot = velocityError;
velocityErrorDot = state.loadAcceleration(:) - desired.acceleration(:);
positionIntegralDot = saturatedIntegralRate(positionIntegral, positionIntegralRate, ...
    cfg.loadController.integralLimit);
desiredForceDot = m0 * (-cfg.loadController.kx .* positionErrorDot ...
    - cfg.loadController.kv .* velocityErrorDot ...
    - cfg.loadController.ki .* positionIntegralDot + desired.jerk(:));

C = R0.' * R0d;
CDot = -hat(Omega0) * C + C * hat(Omega0d);
Q = R0d.' * R0;
QDot = -hat(Omega0d) * Q + Q * hat(Omega0);
eR0Dot = 0.5 * vee(QDot - QDot.');
eOmega0Dot = Omega0Dot - CDot * Omega0d - C * Omega0dDot;
feedforwardRate = C * Omega0d;
feedforwardRateDot = CDot * Omega0d + C * Omega0dDot;
desiredMomentDot = -cfg.loadController.kR .* eR0Dot ...
    - cfg.loadController.kOmega .* eOmega0Dot ...
    + hat(feedforwardRateDot) * J0 * feedforwardRate ...
    + hat(feedforwardRate) * J0 * feedforwardRateDot ...
    + J0 * (CDot * Omega0dDot + C * Omega0dDDot);
if isfield(cfg.loadController, 'yawChannelEnabled') ...
        && ~cfg.loadController.yawChannelEnabled
    desiredMomentDot(3) = 0;
end

% 先把惯性系合力变到负载体系，再求导：
% d(R0' Fd)/dt = R0' Fd_dot - hat(Omega0) R0' Fd。
rhsDot = [R0.' * desiredForceDot - hat(Omega0) * (R0.' * desiredForce); ...
    desiredMomentDot];
muBodyDot = P.' * (PPt \ rhsDot);
muDesiredDot = zeros(3, n);
muDesiredBodyDot = reshape(muBodyDot, 3, n);
for i = 1:n
    bodyIndex = 3 * (i - 1) + (1:3);
    muDesiredDot(:, i) = R0 * (hat(Omega0) * muBody(bodyIndex) ...
        + muDesiredBodyDot(:, i));
end

Omega0DotCmd = cfg.payload.inertia \ ...
    (desiredMoment - cross(Omega0, J0 * Omega0));
Omega0DotCmdDot = J0 \ (desiredMomentDot ...
    - cross(Omega0Dot, J0 * Omega0) - cross(Omega0, J0 * Omega0Dot));

x0ddMinusGE3Dot = desiredForceDot / m0;
forceCommandDot = zeros(3, n);
rhoAll = cfg.payload.attachPoints;
qAll = state.linkUnits;
qDotAll = state.linkRates;
for i = 1:n
    q = normalizeVector(qAll(:, i));
    qDot = qDotAll(:, i);
    rho = rhoAll(:, i);
    muId = muDesiredAll(:, i);
    muIdDot = muDesiredDot(:, i);
    qid = desiredLinkAll(:, i);
    % 虽然主控制器把 omega_id 前馈项关闭，但 e_q = q_id x q_i
    % 本身仍随 q_id(t) 变化，所以 d(e_q)/dt 必须使用解析 qidDot。
    muNorm = norm(muId);
    if muNorm > cfg.loadController.forceNormEpsilon
        qidDot = -(eye(3) - qid * qid.') * muIdDot / muNorm;
    else
        qidDot = zeros(3, 1);
    end

    Omega = linkAngularVelocityAll(:, i);
    ai = accelerationAll(:, i);
    aiDot = x0ddMinusGE3Dot ...
        + R0 * (cross(Omega0, cross(Omega0, cross(Omega0, rho))) ...
        + cross(Omega0Dot, cross(Omega0, rho)) ...
        + cross(Omega0, cross(Omega0Dot, rho))) ...
        - R0 * (hat(Omega0) * hat(rho) * Omega0DotCmd ...
        + hat(rho) * Omega0DotCmdDot);

    q = normalizeVector(q);
    qProjector = q * q.';
    muI = qProjector * muId;
    muIDot = (qDot * q.' + q * qDot.') * muId + qProjector * muIdDot;
    % 当前拍的绳索角加速度由论文式 (7) 直接得到；不能读取上一拍的
    % derivative.linkAngularAccelerations，否则解析链会多出一个 dt 延迟。
    omegaDot = cross(q, ai) / l - cross(q, perpendicularForceAll(:, i)) / (m * l);
    omegaSquared = dot(Omega, Omega);
    omegaSquaredDot = 2 * dot(Omega, omegaDot);

    uParallelDot = muIDot ...
        + m * l * (omegaSquaredDot * q + omegaSquared * qDot) ...
        + m * (qDot * (q.' * ai) ...
        + q * (qDot.' * ai + q.' * aiDot));

    eq = cross(qid, q);
    eqDot = cross(qidDot, q) + cross(qid, qDot);

    % 必须和上面的绳向环使用同一条积分器语义：积分器触及限幅且
    % 误差继续把它推向限幅方向时，实际积分导数为 0；否则为 eq。
    integralRateActual = saturatedIntegralRate(linkIntegralAll(:, i), eq, ...
        cfg.linkController.integralLimit);

    % 主控制器的 USE_QID_DOT=false 使 omega_id=0；因此 e_omega_i=omega_i，
    % 但上面的 e_q 导数仍保留 qidDot。
    correctionDot = -cfg.linkController.kq * eqDot ...
        - cfg.linkController.komega * omegaDot ...
        - cfg.linkController.kqIntegral * integralRateActual;
    qHat = hat(q);
    qHatDot = hat(qDot);
    uPerpDot = m * l * (cross(qDot, correctionAll(:, i)) ...
        + cross(q, correctionDot)) ...
        - m * (qHatDot * qHat + qHat * qHatDot) * ai ...
        - m * qHat * qHat * aiDot;
    forceCommandDot(:, i) = uParallelDot + uPerpDot;
end
end

function rate = saturatedIntegralRate(integralState, rawRate, limit)
%SATURATEDINTEGRALRATE 连续时间积分器的限幅一致导数。
% 离散控制器先执行 I_dot = rawRate，再把 I 截到 [-limit, limit]。
% 在已经触及上限且 rawRate 仍为正（或触及下限且 rawRate 仍为负）时，
% 截幅后的实际导数为 0；其余情况下仍为 rawRate。
if isempty(limit)
    rate = rawRate;
    return
end
limit = limit(:);
if numel(limit) == 1
    limit = repmat(limit, size(integralState));
elseif numel(limit) ~= numel(integralState)
    error('crazyflie_slung_controller:IntegralLimitDimension', ...
        '积分限幅必须是标量或与积分状态同维的向量。');
else
    limit = reshape(limit, size(integralState));
end
if any(~isfinite(limit(:)))
    rate = rawRate;
    return
end
rate = rawRate;
validLimit = limit > 0;
upperActive = validLimit & integralState >= limit & rawRate > 0;
lowerActive = validLimit & integralState <= -limit & rawRate < 0;
rate(upperActive | lowerActive) = 0;
end

% ======================================================================
function OmegaC = analyticOmegaC(Rc, ui, uiDot, b1d, b1dDot, cfg)
%ANALYTICOMEGAC 由 R_c = [b1c b2c b3c] 的**列向量导数**解析求 Omega_c = vee(R_c'RcDot)。
% 与 v2 的 analyticOmegaC **逐行同构**，只有 b3cDot 的种子不同：
%   v2   ：b3cDot = -(I - b3c b3c')*A_dot/||A||，其中 A_dot 可由位置环解析得到；
%   本工程：u_i^cmd 的导数由 analyticForceCommandDerivative 按
%          Fd/Md -> mu_id -> q_id -> a_i -> u_parallel/u_perp 的链式法则给出。
b3c = Rc(:, 3);
forceNorm = norm(ui);
if forceNorm < cfg.loadController.forceNormEpsilon
    OmegaC = zeros(3, 1);
    return
end
b3cDot = -(eye(3) - b3c * b3c.') * uiDot / forceNorm;

% b1c 是 b1d 在垂直于 b3c 平面上的归一化投影（与 (37) 一致）
b1c = Rc(:, 1);
projection = (eye(3) - b3c * b3c.') * b1d;
if norm(projection) < cfg.loadController.forceNormEpsilon
    % 航向投影退化时控制器用了备用方向；此处令其导数为零，
    % 避免在奇异点附近产生异常大的角速度前馈（与 v2 一致的处理）。
    b1cDot = zeros(3, 1);
else
    projectionDot = -(b3cDot * b3c.' + b3c * b3cDot.') * b1d ...
        + (eye(3) - b3c * b3c.') * b1dDot;
    b1cDot = (eye(3) - b1c * b1c.') * projectionDot / norm(projection);
end

% b2c = b3c x b1c 的导数
b2c = Rc(:, 2);
b2cDot = cross(b3cDot, b1c) + cross(b3c, b1cDot);
RcDot = [b1cDot, b2cDot, b3cDot];
OmegaHat = Rc.' * RcDot;
OmegaHat = 0.5 * (OmegaHat - OmegaHat.');
OmegaC = vee(OmegaHat);
end

% ======================================================================
function v = so3Log(R)
%SO3LOG SO(3) 对数映射：旋转矩阵 -> 旋转向量（角度 = 向量模长）。与 v2 一致。
R = projectSO3(R);
cosTheta = min(max((trace(R) - 1) / 2, -1), 1);
theta = acos(cosTheta);
if theta < 1e-7
    v = 0.5 * vee(R - R.');
elseif abs(pi - theta) < 1e-5
    % 接近 180 度：用 (R + I)/2 的最大特征向量作为旋转轴
    [V, D] = eig((R + eye(3)) / 2);
    [~, idx] = max(real(diag(D)));
    axisVector = real(V(:, idx));
    axisVector = axisVector / max(norm(axisVector), eps);
    v = theta * axisVector;
else
    v = theta / (2 * sin(theta)) * vee(R - R.');
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
%OUTWARDINTERNALBIAS 构造不改变总 wrench 的径向内部张力。
n = size(rhoAll, 2);
desired = zeros(3, n);
for i = 1:n
    radial = [rhoAll(1, i); rhoAll(2, i); 0];
    radialNorm = norm(radial);
    if radialNorm < 1e-12
        angle = 2 * pi * (i - 1) / max(n, 1);
        radial = [cos(angle); sin(angle); 0];
    else
        radial = radial / radialNorm;
    end
    desired(:, i) = fraction * mass * gravity / sqrt(n) * radial;
end

% 投影到 null(P)，所以 P*bias = 0，合力和合力矩严格不变。
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
