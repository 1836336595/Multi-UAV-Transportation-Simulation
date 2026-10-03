function cfg = crazyflie_slung_parameters(userCfg)
%CRAZYFLIE_SLUNG_PARAMETERS 多机（n 架）协同吊运刚体负载仿真的集中参数文件。
%
% 直接对应 Lee 2014 / Lee 2018 的论文设定：n 架四旋翼各自通过一根无质量**绳索**
% 吊在同一刚体负载的不同挂点上，协同完成负载的位姿控制。
% ★ 论文原文用词是 "massless links"、标题是 "Cable-Suspended Rigid Body"，
%   其 Remark 1 明确要求"每根缆的张力为正"——所以物理上是**绳（cable）**，
%   不是刚性连杆。绷紧的绳与无质量刚性连杆在力学上完全等价（都只沿 q_i 传轴向力），
%   差别只有一条：绳多一个单边约束 mu_i >= 0（只能拉，不能推）。
%   TAUT_ACTIVE 段要求 mu_i > 0；起飞/降落的 SLACK 段按零张力处理，详见 README。
%   负载：    位置 x0、姿态 R0、质量 m0、惯量 J0、挂点 rho_i (i = 1..n)
%   第 i 机： 质量 m_i、惯量 J_i、绳长 l_i、绳向 q_i（由该机指向负载）
%   几何约束：x_i = x0 + R0 rho_i - l_i q_i
%   控制输入：{f_i, M_i}，本实现把它改写为 {推力百分比, 机体角速度指令}
%             （对应 cflib 的 send_setpoint_manual(rate=True)）
%
% 用法：
%   cfg = crazyflie_slung_parameters();
%   cfg = crazyflie_slung_parameters(userCfg);   % userCfg 覆盖默认值
%
% 坐标系与符号约定（完全沿用论文）：
%   * 惯性系第三轴 e3 = [0;0;1] 指向重力方向，因此 z 轴向下为正；
%   * R0, R_i 均为"机体系 -> 惯性系"的旋转矩阵，R_dot = R * hat(Omega)；
%   * 负载位置 x0，绳索单位向量 q_i（由第 i 架四旋翼指向负载）；
%   * 四旋翼推力为 -f_i R_i e3（f_i > 0 时指向机体上方），故负载吊在下方。

if nargin < 1 || isempty(userCfg)
    userCfg = struct();
end

% ---------------------------------------------------------------- 仿真设置
% ★ 本分支是**定高**工况（无水平机动），默认时长 30 s。
%   历史说明：静态工况下不存在真稳态，位置误差在 t = 30 s 取最小 3.63 mm 后
%   回升、77.2 s 完全发散（旧版本 yaw 符号和反馈配置错误），故 30 s 作为
%   默认展示窗口；当前 yaw 通道默认开启并单独记录跟踪误差。
%   dt = 0.002 s 对应 500 Hz。
cfg.simulation = struct(...
    'duration', 30.0, ...          % 定高悬停 30 s（悬停收敛需要更长时间）
    'dt', 0.002, ...               % 控制/积分步长 [s]
    'maxPendulumRate', 12.0, ...   % 绳索角速度安全上限 [rad/s]，仅数值保护
    'maxBodyRate', 12.0);          % 机体角速度安全上限 [rad/s]，仅数值保护

% ------------------------------------------------------------------ 负载
% 悬挂一块长方体平台。n 架四旋翼分别挂在平台上三个不共线的挂点上。
cfg.payload = struct();
cfg.payload.mass = 0.080;                          % m0，负载质量 [kg]
cfg.payload.size = [0.08; 0.06; 0.05];            % [长; 宽; 厚] 长方体外形 [m]
% 几何派生量在 mergeStruct 后按 size 重新计算。mass 默认独立于尺寸；
% 若要模拟同材料物品随体积变重，可在 userCfg.payload 中显式提供 density，
% 并省略 mass，此时质量按 density * volume 自动计算。
% 下面两个字段是**无量纲几何模板**，实际挂点由 payload.size 自动缩放：
%   x = +/- a/2，y = 0 或 +/- b/2，z = -c/2（上表面，z 轴向下）。
%   修改 payload.size 后，挂点仍然落在对应的物品棱/顶面边界上。
cfg.payload.attachFractions = [ 0.5, -0.5, -0.5; ...
                                0.0,  0.5, -0.5; ...
                               -0.5, -0.5, -0.5];
% 惯量按均质长方体自动计算：
%   Ixx = m(b^2+c^2)/12, Iyy = m(a^2+c^2)/12, Izz = m(a^2+b^2)/12
cfg.payload.inertia = boxInertia(cfg.payload.mass, cfg.payload.size); % J0 [kg*m^2]

% ★★★ 负载**转动阻尼**（物理项，非控制项）—— 2026-09-21 新增，必须保留
%   背景：用户要求的挂点几何（一边中点 + 对边两顶点）可能让挂点几何中心偏离负载质心，
%         偏移量由 payload.size 和 attachFractions 自动决定。该偏移带来一个
%         "负载偏航 <-> 绳索扭转"的耦合模态，
%         频率 ω ≈ sqrt(m0 g offset^2 / (L J0z)) ≈ 2.2 rad/s。
%         而负载偏航通道必须通过绳索倾斜来交付。修正动力学中的
%         hat(rho_i) 力矩符号后，yaw 反馈可以恢复为负反馈；转动阻尼仍保留，
%         用于耗散离散执行器和挂点偏心引起的残余振荡。
%   解决：把**真实存在的转动阻尼**补进模型。它是**外部力矩**，直接加在负载转动
%         方程右端，**不经控制器、不进张力分配**，因此不会触发上述正反馈。
%   取值依据（物理量级估算，不是凑参）：参考负载尺寸由
%   cfg.payload.rotationalDampingReferenceSize 给出，阻尼随后按水平面积平方缩放。
%         气动阻尼力矩约 (1/8) rho Cd a^4 omega^2；取 rho=1.225、Cd≈1.1、
%         ω=3 rad/s 得 ≈ 2.4e-3 N*m，等效线性系数 c ≈ 8e-4 N*m*s/rad。
%         取 1.0e-3（同一量级，略偏保守）。
%   实测效果（镜像）：max|Omega0_z| 由 4.058 降到 0.19（约 21 倍），
%         位置跟踪峰值 705.0 mm **完全不变**，推力峰值 69.1% 不变。
%         c 的有效区间很宽：5e-4 ~ 1e-2 都稳定；>= 5e-2 才会发散（过阻尼+离散化）。
cfg.payload.rotationalDamping = 1.0e-3;   % [N*m*s/rad] 负载转动阻尼系数
cfg.payload.rotationalDampingReferenceSize = cfg.payload.size;
% 挂点（负载体系坐标，z 向下为正，故 z=-c/2 是上表面）。列 = 挂点，
% 即论文的 rho_i。实际值在文件末尾根据 payload.size 生成。
%
% ★ 挂点按用户要求放在**物品边缘**（不再放物品中间）：三个点构成一个等腰直角三角形，
%   其中一点是 x=+a/2 边的中点，另两点是 x=-a/2 边的两个顶点：
%     rho_1 = ( +a/2,  0,   -c/2 )
%     rho_2 = ( -a/2, +b/2, -c/2 )
%     rho_3 = ( -a/2, -b/2, -c/2 )
%   其中 [a;b;c] = payload.size，因此修改长宽高后挂点自动保持在边界上。
%
% ★★ 这个几何的两个后果（已数值核算，见 README §11）：
%   1) **分配矩阵条件数变好**：cond(P*P') 从 276.8 降到 150.2（挂点三角形更大，
%      力臂更长）。rank(P) 仍为 6，可解性不受影响。
%   2) **出现静载不对称 2:1:1**：挂点质心在 (-0.0333, 0)，相对负载质心偏了 0.0333 m，
%      于是悬停时三根绳必须给出**不等张力**才能平衡力矩：
%         T = [0.3924, 0.1962, 0.1962] N   （旧的内接等边三角形是 0.2616 N 三等分）
%      对应推力占比 [53.3%, 38.6%, 38.6%]（旧值 [43.5% ×3]）。
%      ⇒ 机 1 的推力需求提高，机 2/3 降低；最小张力 0.1962 N（旧悬停值的 75%），
%        绳索绷紧裕度相应减小，**必须靠自检确认全程无松弛**。
cfg.payload.attachPoints = bsxfun(@times, cfg.payload.size(:), ...
    cfg.payload.attachFractions);

% ------------------------------------- Crazyflie 2.1 Brushless
% ★★ 四旋翼**不做刚体姿态建模**（与 v2 的 crazyflie_ctbr_simulation.m 一致）：
%   机体转动惯量既不需要、也不影响任何结果。原因可以自己验算 ——
%   速率环给出的角加速度是 α_cmd，若再把它包成力矩
%        M_i = J_i·α_cmd + Ω_i × J_i Ω_i
%   而姿态方程又是
%        Ω̇_i = J_i⁻¹ ( M_i − Ω_i × J_i Ω_i )
%   代回去恰好得到 **Ω̇_i = α_cmd**，J_i 与叉乘项完全抵消。
%   ⇒ 直接用"角速度指令 → 一阶速率环 → 角速度 → 积分得到姿态"：
%        Ω̇_i = rateLoop.bandwidth .* (Ω_cmd_i − Ω_i)
%        R_i ← R_i · exp(hat(Ω_i)·dt)
%   这也是真实 Crazyflie 的用法：飞控只接受**角速度指令**
%   （cflib 的 send_setpoint(..., rate=True)），力矩由固件自己的速率环产生，
%   外部无需知道机体惯量。
%
%   质量：含电池、无桨叶保护罩的标称起飞质量 0.0325 kg（Bitcraze 官方规格）。
%   最大推力：4 x 34 gf = 136 gf ≈ 1.334 N（Bitcraze 官方参数），只用作推力饱和限幅。
%   推力时间常数：一阶执行器响应；真实飞行前应通过台架辨识替换。
%   臂长 / 旋翼半径 / 机体外形 / 视觉转速：**只给三维显示与碰撞半径用**，不进动力学。
cfg.vehicle = struct(...
    'count', 3, ...                                % n，四旋翼数量（须 >= 3，见文件末尾检查）
    'mass', 0.0325, ...                            % m_i [kg]
    'gravity', 9.81, ...                           % g [m/s^2]
    'armLength', 0.0465, ...                       % 电机轴到质心距离 [m]（显示 + 碰撞半径）
    'maxTotalThrust', 1.3344, ...                  % 单机 4 电机最大总推力 [N] = 136 gf（推力饱和限幅）
    'thrustTimeConstant', 0.012, ...               % 推力一阶响应时间常数 [s]
    'bodySize', [0.050; 0.050; 0.014], ...         % 机体中心板外形 [m]（仅三维显示）
    'rotorRadius', 0.023, ...                      % 单个旋翼半径 [m]（显示 + 碰撞半径）
    'rotorSpinHz', 4.0);                           % 动画里旋翼视觉转速 [Hz]（仅显示）

% ------------------------------------------------------ 绳索（只受拉、无质量、绷紧时与刚性连杆力学等价）
% ★ 物理是"绳"：只能受拉（单边约束 mu_i >= 0），不能受压。
%   绷紧时它的力学与无质量刚性连杆**完全相同**，所以论文的动力学与控制律无需改动；
%   只有当 mu_i 下探到 0 时绳才松弛、负载自由落体，那时才需要额外处理。
%   TAUT_ACTIVE 段的张力必须保持正值；地面起飞前不调用该绷紧段模型。
cfg.link = struct(...
    'length', 0.65, ...                            % l_i [m]
    'count', 3, ...                                 % 与 cfg.vehicle.count 保持一致
    'allowTiltedCables', true, ...                 % 允许 q_i 不平行于 e3
    'initialOutwardOffset', NaN, ...                % 初始无人机相对挂点的外张距离 [m]
    'vehicleClearance', 0.05, ...                    % 机体中心额外安全间隙 [m]
    'minInFlightTiltRatio', 0.20);                 % 空中绳向最小倾角比例（相对绳长）

% --------------------------------------------------------- 起飞/降落混合阶段
% 这组参数只控制仿真中的地面接触、松弛绳和绷紧过渡，不改变 Lee 的绷紧段
% 动力学。没有绳端拉力传感器时，TAKEUP -> TAUT_RAMP 使用几何距离判据：
%   d_i = ||(x_0 + R_0 rho_i) - x_i||,
% 再配合 epsilonOn/epsilonOff 迟滞，避免定位噪声造成状态抖动。
% z 轴沿用论文约定（向下为正），所以 groundZ=0、负载中心在地面时为
% groundZ - payload.size(3)/2，显示时再取反为正高度。
cfg.takeoff = struct(...
    'enabled', true, ...
    'landingEnabled', true, ...
    'groundZ', 0.0, ...                         % 地面 z 坐标（惯性系）
    'vehicleGroundClearance', 0.045, ...        % 机体中心离地高度 [m]
    'groundFrictionMu', 0.50, ...               % 负载触地时的库仑摩擦系数
    'groundRotationalBrake', 40.0, ...          % 触地时角减速度上限 [rad/s^2]
    'groundContactTolerance', 0.002, ...        % 判定"仍在地面接触"的竖直余量 [m]
    'groundRadialOffset', NaN, ...              % 起降阶段相对负载中心的安全外张 [m]
    'independentHoverHeight', 0.30, ...         % 起飞段目标离地高度 [m]
    'takeoffDuration', 2.5, ...                % 独立起飞到收紧高度 [s]
    'takeupDuration', 2.5, ...                  % 从收紧高度缓慢接近绳长 [s]
    'preTensionSlack', 0.0, ...                 % 收紧末端保留的绳长余量 [m]（★ 必须 0，见下）
    'takeupSnapWarn', 0.010, ...                % 交接瞬移超过它就告警 [m]
    'epsilonOn', 0.030, ...                     % 进入绷紧候选的距离余量 [m]
    'epsilonOff', 0.060, ...                    % 释放判据的距离余量 [m]
    'confirmTime', 0.30, ...                    % 距离条件持续时间 [s]
    'tensionRampTime', 1.50, ...                % 张力软建立时间 [s]
    'referenceLiftTime', 3.00, ...              % 张力建好后参考抬升到目标的时间 [s]
    'releaseSnapTolerance', 0.005, ...          % 准许"落地点释放"的离地余量 [m]
    'landingDuration', 4.0, ...                 % 末段受控下降 + 独立降落 [s]
    'landingApproachFraction', 0.55, ...        % 前一部分仍由绷紧动力学下降
    'independentPositionKp', [3.0; 3.0; 4.0], ...
    'independentPositionKv', [3.12; 3.12; 3.60], ...
    'independentIntegralGain', [0.80; 0.80; 0.80], ...
    'independentMaxFeedbackAcceleration', [6.0; 6.0; 8.0], ...
    'independentIntegralLimit', [0.20; 0.20; 0.20], ...
    'independentIntegralGate', 0.050, ...        % 位置误差小于它才积分（抗饱和）[m]
    'independentMaxBodyRate', [5.0; 5.0; 3.5], ...
    'independentHeading', [1; 0; 0]);

% ★★ 为什么 takeoff.preTensionSlack 必须是 0（2026-09-29）★★
%   本模型**只能表示绷紧的绳**（松弛段按设计未实现）。而"收紧末端保留 σ 的绳长余量"
%   在数学上等价于"进入绷紧段时把 σ 一次性收掉" —— 实测 σ = 12 mm 时表现为
%   **三架无人机在 2 ms 内同时瞬移 17 mm**（12 mm 余量 + 5 mm 独立控制器稳态余差，
%   等效速率 8.5 m/s），直接激励绳向环。
%   ⇒ 正确的交接点是"绳刚好拉直、张力为零"，即 σ = 0。
%   σ > 0 仍可设置（校验允许 [0, L)），但会在交接瞬间产生同样大小的瞬移，
%   此时 simulation.m 会按 takeoff.takeupSnapWarn 发出 warning。
%   残余量 ≈ 独立控制器的稳态余差（实测 ~5 mm）；要再压小需提高
%   independentIntegralGain（积分收敛时间 ≈ 2*Kp/Ki，当前 0.8 ⇒ 约 10 s，
%   比 TAKEUP 的 2.4 s 长，所以收敛不完）。

% ------------------------------------------- 负载位置/姿态外环（论文 (20)-(21)）
% ★ 重要：论文 (20) 式的等效质量是 **m0**（负载质量），不是 (m0 + sum m_i)。
%   这一点容易写错。验证方式：代入 (17) 与 (5)/(6) 化简可得 (18)
%        m0 (x0_ddot - g e3) = sum_i mu_i
%   即虚拟控制（绳索张力）之和恰好平衡负载自身的重力。
%   物理直觉：绳索张力只需要托住负载 m0 g；四旋翼自身的重量由各自的推力
%   通过 (17) 式中的 m_i q_i q_i' a_i 项补偿，不进入 mu。
%   若误把 Mq = (m0 + sum m_i) I + sum m_i q_i q_i' 用到这里，悬停推力会被
%   多算 sum(m_i) g，负载会持续下沉直到积分器补上，出现明显稳态高度偏差。
%
% 增益选取依据（内环带宽必须显著高于外环）：
%   速率环 35 rad/s  >  摆频 sqrt(g/l) = 5.29 rad/s  >  负载回路 1.73 rad/s
%
% ★ ki 是稳态精度的主控，且**加大 ki 在三机情形下没有副作用**（_sweep_3drone.py）：
%     ki = 0.60 -> 稳态位置 23.99 mm，绳向误差 1.410 deg，推力峰值 87.3 %
%     ki = 0.90 -> 稳态位置 21.28 mm，绳向误差 1.411 deg，推力峰值 87.3 %
%     ki = 1.20 -> 稳态位置 17.72 mm，绳向误差 1.411 deg，推力峰值 87.3 %
%     ki = 1.60 -> 稳态位置 13.56 mm，绳向误差 1.412 deg，推力峰值 87.3 %
%   位置误差单调下降而绳向误差与推力峰值**完全不变**，说明 ki 只作用于
%   负载平动通道，不与摆环/姿态环争夺权限。故取 1.60。
%
% ★★ 负载姿态增益必须按 J0 缩放（这是本仿真最容易踩的坑之一）：
%   论文 (21) 的 Md = -kR e_R0 - kOmega e_Omega0 是"未除惯量"的力矩形式，
%   闭环二阶特性为
%       J0 * Omega0_dot = -kR e_R0 - kOmega e_Omega0
%   即   omega_n = sqrt(kR / J0)，  zeta = kOmega / (2 sqrt(kR J0))
%   本负载 J0 ≈ 2.69e-4 kg*m^2 极小。若沿用"按 J0 = 1 调出"的增益
%   （kR = 0.55），omega_n = sqrt(0.55/2.69e-4) = 45 rad/s ≈ 7.2 Hz 尚可，
%   但阻尼项 kOmega = 0.35 给出的等效带宽 kOmega/J0 = 1300 rad/s ≈ 207 Hz，
%   已逼近 500 Hz 采样率的奈奎斯特边界，姿态环数值发散。
%   下面的取值按目标带宽与阻尼比反解：
%       平移通道：omega_n = 2*pi*6.0 = 37.7 rad/s，zeta = 0.9
%     kR     = J0 .* omega_n^2
%     kOmega = 2 * zeta * omega_n .* J0
%
% ★ 负载偏航（第 3）通道默认开启，但采用远低于 roll/pitch 的带宽。
%   论文中的 yaw 控制力矩 Md(3) 由分配矩阵 P 分解到各根绳索，再由绳向环
%   产生所需的水平张力分量。n = 3 且挂点不共线时 rank(P) = 6，理论上
%   合力和三轴合力矩（包括 yaw）均可分配。yaw 交付能力取决于挂点力臂、
%   绳索倾角和执行器裕度，因此带宽不能直接照搬 roll/pitch。
%   代码曾把式（7）中的 hat(rho_i) 误写成 hat(rho_i)'，该转置会将绳索
%   力矩整体反号，使 yaw 反馈变成正反馈；现已修正为论文原符号。
payloadInertia = cfg.payload.inertia;
wnLoad = 2 * pi * [6.0; 6.0; 0.45];      % yaw 低带宽，避免压过绳索摆动模态
zetaLoad = 0.90;                        % 目标阻尼比
kRLoad = diag(payloadInertia) .* wnLoad.^2;
kOmegaLoad = 2 * zetaLoad * wnLoad .* diag(payloadInertia);

% ★ 手动指定负载姿态增益（留空 [] = 用下面的自动计算；
%   放在 struct 外面是为了不打断续行，也便于静态检查识别字段）
cfg.loadController = struct(...
    'kx', [3.0; 3.0; 3.75], ...       % 定高增益（omega_n ≈ 1.73~1.94 rad/s）
    'kv', [3.12; 3.12; 3.12], ...     % 定高增益（zeta ≈ 0.90）
    'ki', [1.60; 1.60; 1.60], ...     % 定高增益（受限，抗常值扰动）
    'kR', kRLoad, ...                 % 负载姿态增益 [N*m/rad] = J0 .* wn^2
    'kOmega', kOmegaLoad, ...         % 负载角速度增益 [N*m*s/rad] = 2 zeta wn J0
    'c1', 0.50, ...                   % 积分器交叉项系数（抑制饱和下的过冲）
    'integralLimit', [0.50; 0.50; 0.50], ...   % 积分饱和限幅
    'forceNormEpsilon', 1e-9, ...      % ||Fd|| 数值保护阈值
    'manualKR', [], ...               % 手动 kR     [N*m/rad]  空=自动
    'manualKOmega', []);              % 手动 kOmega [N*m*s/rad] 空=自动
% ★★★ 位置增益是**绕八字工况专用**的重标定值，与静态悬停工况不同 ★★★
%   静态悬停时用 [3.0; 3.0; 3.75] 就够（误差收敛到 3.63 mm）。
%   但绕八字时参考点在持续运动，位置环必须"追得上"，否则留下
%   相位滞后型的跟随误差。实测（_diag_fe_tune4.py / _diag_fe_tune5.py，A=±1.5 m）：
%
%     KX     wn=sqrt(KX)   相对绳摆频率        环绕段误差均值   末点残差   稳定性
%     3.0    1.73 rad/s    0.33 x            287.0 mm        104.7 mm  稳定
%     5.0    2.24 rad/s    0.42 x            216.9 mm        130.9 mm  稳定
%     8.0    2.83 rad/s    0.53 x            162.7 mm        103.3 mm  稳定
%    10.0    3.16 rad/s    0.60 x            142.8 mm         96.4 mm  稳定
%    14.0    3.74 rad/s    0.71 x            123.5 mm         78.7 mm  稳定 ← 采用
%    20.0    4.47 rad/s    0.84 x            ——               ——        11.10 s 发散
%    24.0    4.90 rad/s    0.93 x            ——               ——         3.83 s 发散
%    28.0    5.29 rad/s    1.00 x（匹配）    ——               ——         2.42 s 发散
%
%   ★ 一个反直觉的结论：**让位置环带宽"匹配"绳摆频率反而会发散**。
%     纯直觉会说"带宽匹配才有最好的能量传递"，但这里的耦合是**寄生**的：
%     位置环越快，它越会把负载的水平运动"硬"地转成绳索摆动，而绳向环
%     （kq = 55）跟不上这个加快的激励，摆角迅速积累到翻转发散。
%     所以最优带宽是**绳摆频率的 0.6~0.7 倍**（约 -3 dB 处仍能跟随，但
%     不与摆模态共振），而不是 1.0 倍。
%   ★ 边界很陡：14 -> 20 之间只有一步之差就发散。因此 KX = 14 是
%     "贴着边界取"的最优值，靠下面三项裕度守住：
%       a. 参考剖面用 5 次多项式（两端速度/加速度为 0，消除参考加速度尖峰）
%       b. USE_QID_DOT = false（去掉 q_id_dot 差分噪声，见 controller.m 说明）
%       c. KV/KX = 0.62（zeta ≈ 0.90，给摆模态留阻尼裕度）
%   ★ KV 的取法：闭环二阶特性 wn = sqrt(KX)，zeta = KV / (2 sqrt(KX))。
%     取 zeta = 0.90 -> KV = 1.8 * sqrt(KX) = 1.8 * 3.742 = 6.73。
%     实测 0.62 * KX = 8.68 更好（等效 zeta = 1.16，略过阻尼），
%     因为过阻尼能进一步压制摆模态的相位滞后。24c 的阻尼细扫确认
%     KV/KX 在 0.62 附近最稳，偏离到 0.50 或 0.90 都提前发散。
%   ★ ki 同步从 1.60 提到 3.00：KX 提高后同样的 ki 显得偏软，
%     末点残差从 123 mm 降到 78.7 mm。ki >= 6.0 会让积分在环绕段
%     的往复激励下发散（_diag_fe_tune3.py 实测 KI=6.0 于 11.30 s 发散）。
%   竖直通道 kx(3) = 1.25 * 14.0 = 17.5，与水平通道同一比例（原来 3.75/3.0 = 1.25）。
% 记录设计意图，便于复现与调参
% yaw 通道开启；其带宽远低于 roll/pitch，避免用高增益直接激励绳向环。
cfg.loadController.designBandwidthHz = wnLoad / (2 * pi);
cfg.loadController.yawChannelEnabled = true;
cfg.loadController.designDampingRatio = zetaLoad;

% 姿态力矩『张力预算』上限（2026-09-28 新增）
%   病根：kR = J0*wn^2 使期望力矩 Md 正比于 J0（尺寸的平方）；
%   而把它交付到负载上靠各绳的差分张力 delta_mu ≈ Md/(n*rho)。
%   可用张力是 m0*g/n（常数）⇒ delta_mu/可用 正比于尺寸。
%   实测（mass 固定、尺寸放大 k 倍）：k=1 → 14.7%、k=2 → 29.5%、
%   k=2.5 → 36.8%、k=4 → 58.9% ⇒ 越大越接近吃光张力预算
%   ⇒ 绳趋松弛、绳向环追不上 ⇒ 发散。这是改尺寸后不稳的主因。
%   修法：由张力预算反推 kR 上限，大负载自动降带宽换可行性。
cfg.loadController.attitudeMomentBudget = 0.35;   % 允许力矩占 m0*g*rhoTyp 的比例
cfg.loadController.attitudeMomentRefError = 0.10; % 折算用参考姿态误差 [rad]

% --------------------------------------------- 张力分配（论文 (13)(22)-(25)）
% 需要的合力/力矩 [Fd; Md] 通过分配矩阵 P 的伪逆分配到 n 根绳索的张力上：
%     P = [ I, I, ..., I ; hat(rho_1), hat(rho_2), ..., hat(rho_n) ]   (6 x 3n)
%     [mu_1d; ...; mu_nd] = diag[R0,...,R0] * P' * (P*P')^-1 * [R0'*Fd; Md]
% ★ 可解性条件：rank(P) = 6。P 为 6 x 3n，故必须满足
%     * n >= 2（3n >= 6），但 n = 2 时存在零空间
%        spanned by [(rho_1-rho_2); (rho_2-rho_1)]，
%        因此论文 (2014) 明确要求 **n >= 3** 才能保证张力分配对
%        所有 [Fd; Md] 唯一可解、且不产生内力环流。
%     * n 个挂点不能共线（否则 rank[hat(rho_1) ... hat(rho_n)] < 3）。
%   单架四旋翼（n = 1）时 3n = 3 < 6，无法同时满足 3 个力分量 + 3 个力矩分量，
%   因此论文的多机框架在 n = 1 下退化。
%   本仿真的实测值：n = 3 时 cond(P*P') = 276.87，rank(P) = 6，条件良好。
cfg.allocation = struct(...
    'pinvTolerance', 1e-9, ...        % P*P' 求逆前的条件数检查阈值
    'checkRank', true, ...             % 是否在启动时检查 rank(P) 并给出提示
    'outwardBiasFraction', 0.20, ...  % 内部张力外张偏置 / (m0*g)
    'outwardBiasMax', 0.12);           % 每根绳外张偏置的 RMS 上限 [N]

% ---------------------------------------- 绳向环（论文 (26)-(28)）
% 论文 (27)：u_perp_i = m_i l_i hat(q_i){ -kq e_qi - kw e_wi - (q_i.omega_id) q_id_dot
%                                        - hat(q_i)^2 omega_id_dot } - m_i hat(q_i)^2 a_i
% 摆的自然频率 sqrt(g/l) = 5.294 rad/s，摆回路带宽需高于它：
% 取 kq = 55 -> omega_n = 7.42 rad/s = 1.40 倍自然频率（否则摆角会翻转发散）。
% komega = 20 -> zeta = komega/(2*omega_n) = 1.35（略过阻尼，抑制残余振荡）。
%
% ★★ kq 已处于 3 机情形的稳定上限，**不要再往上加**（_sweep_3drone.py 实测）：
%      kq = 55  -> 收敛（本实现采用）
%      kq = 70  -> 发散于  9.2 s
%      kq = 85  -> 发散于  5.3 s
%      kq = 100 -> 发散于  3.5 s
%      kq = 120 -> 发散于  3.0 s
%   机理：负载偏航角会通过挂点方位变化影响绳向环，
%   挂点方位随之转动，因此**期望绳向本身在低频变化**。提高 kq 会让
%   绳向环更"硬"地追赶这个低频漂移，反而把漂移能量放大进姿态环。
%   这与单机版"kq 越大越好"的直觉相反 —— 多机情形下摆环与负载姿态环
%   通过挂点位置强耦合，不能独立整定。
%   若需要更硬的绳向，应先保证 yaw 参考变化平滑并重新检查张力/推力裕度。
cfg.linkController = struct(...
    'kq', 55.0, ...                   % 摆方向增益 [1/s^2]，omega_n ≈ 7.42 rad/s，★ 稳定上限
    'komega', 20.0, ...               % 摆角速度增益 [1/s]，zeta ≈ 1.35
    'kqIntegral', 0.0, ...            % 摆方向积分增益（默认关闭，见下）
    'integralLimit', [0.30; 0.30; 0.30], ...   % 积分饱和限幅
    'maxRelativeAngle', 1.8);         % e_qi 与 e_wi 的数值保护限幅 [rad]
% ★ 摆环积分 kqIntegral 默认关闭，且**实测证明开了也没用**（_sweep_3drone.py）：
%     kqIntegral = 0    -> 稳态位置 23.99 mm，绳向误差 1.410 deg
%     kqIntegral = 0.3  -> 稳态位置 24.00 mm，绳向误差 1.412 deg
%     kqIntegral = 1.5  -> 稳态位置 24.02 mm，绳向误差 1.420 deg
%   开了积分反而**略微变差**。原因是绳向的残差不是常值扰动（积分器能消），
%   而是由偏航漂移持续注入的低频误差，积分器追不上，只引入相位滞后。

% ------------------------------------------------ 姿态外环（替代原文力矩环）
% 原文 (39)-(40) 直接给出力矩 M_i；本方案改为给出角速度指令 omega_cmd_i，
% 由 Crazyflie 固件速率环（本仿真用一阶闭环模型等效）跟踪。
% 姿态环等效时间常数 tau ≈ (1 + kOmega) / kR。
%
% ★★ kR 从 [30;30;15] 提高到 [240;240;120]（×8）—— 挂点改成"边中点+对边顶点"后必需 ★★
%   新挂点几何使挂点质心偏离负载质心 0.0333 m，悬停时三根绳必须给出
%   **不等张力 2:1:1 = [0.3924, 0.1962, 0.1962] N** 才能平衡力矩。
%   每架机受到的绳索反作用力因此不再对称，机体必须更"硬"地保持姿态，
%   否则推力方向偏离 -> 与绳摆/负载姿态环形成正反馈。
%
%   实测（_diag_attach_retune.py，静态悬停 20 s，仅改 kR 倍率）：
%       ×1  (30,30,15)     悬停 1.202 s 发散 / 八字 1.646 s 发散
%       ×2  (60,60,30)     悬停稳定 15.25 mm / 八字 7.290 s 发散
%       ×3  (90,90,45)     悬停稳定 15.11 mm / 八字 10.200 s 发散
%       ×5  (150,150,75)   悬停稳定 14.69 mm / 八字**稳定** 692.3 mm
%       ×8  (240,240,120)  悬停稳定 14.65 mm / 八字稳定 693.5 mm   ← 采用
%       ×20 (600,600,300)  悬停稳定 14.67 mm / 八字稳定 695.5 mm
%   ★ 悬崖在 ×3 与 ×5 之间；×5~×20 是性能几乎不变的平台区，取 ×8 居中留裕度。
%   ★ 这不是"用增益压住问题"：绳索反作用不对称是真实的物理扰动，
%     提高机体姿态环带宽是对应的、正确的处置。
%
% ★★ 外层姿态环的**结构**（2026-10-03 改为与 v2 同构）：
%     Omega_cmd = R'R_c·Omega_ic            <- 期望姿态角速度**前馈**
%                 - kR  .* e_R
%                 [- kOmega .* (Omega - R'R_c·Omega_ic)]   <- 仅 useRateDamping=true
%                 - kIR .* ∫e_R             <- 外层积分，默认关闭（与 v2 的
%                                              useIntegral=false 一致）
%   ★★ `useRateDamping` 决定要不要那一项**速率阻尼**：
%      false（默认）= **与 v2 完全一致的式子**（v2 没有这项，阻尼交给内环 PI）
%      true         = 论文原式（Lee 2018 (30) 的 e_Omega 项），保留 kOmega
%   `omegaCMethod` 决定前馈 Omega_ic 的计算方式：
%      'analytic'（默认）= 按 Fd/Md -> mu_id -> q_id -> a_i -> u_i^cmd 的完整链式法则求导
%      'filtered_log_difference' = 相邻期望姿态的 SO(3) 对数差分
%      'command_filter' = 先用二阶命令滤波器平滑 u_i^cmd，同时直接得到其导数，
%                         滤波后的力同时用于姿态和推力，保证两者一致
%      'high_gain_observer' = 用二阶高增益非线性微分器估计 u̇_i^cmd
%      'none' = 关闭前馈（Omega_ic ≡ 0）
%   ★ 为什么仍要限幅/低通：前馈的种子（u_i^cmd 的变化）含张力分配/绳向环的高频
%     抖动，裸的 1/dt 差分会把抖动放大成每秒几百弧度的假前馈（1/dt = 500）。
%     本项目实测过该自激（姿态误差 0.6 deg 随时间涨到 8.9 deg、机体角速度
%     长期 4.2 rad/s），所以历史上是把前馈整项置零的。现在按 v2 结构启用，
%     但用 `feedforwardFilterTime` 低通、并用 `feedforwardMaxRate` 限幅。
%     若仍见姿态误差缓慢增长 ⇒ 调大 tau，或设 omegaCMethod = 'none'。
cfg.attitudeController = struct(...
    'kR', [240.0; 240.0; 120.0], ...  % 姿态误差 -> 角速度 的增益 [rad/s]
    'kOmega', [4.0; 4.0; 4.0], ...    % 角速度误差 -> 角速度 的增益（仅 useRateDamping=true 时生效）
    'useRateDamping', false, ...      % 是否加 -kOmega·e_Omega（false = 与 v2 完全一致）
    'omegaCMethod', 'high_gain_observer', ...   % 'analytic' | 'command_filter' | 'high_gain_observer' | 'filtered_log_difference' | 'none'
    'feedforwardFilterTime', 0.020, ... % 前馈的一阶低通时间常数 [s]（0 = 不滤波）
    'feedforwardMaxRate', 20.0, ...   % 前馈原始值限幅 [rad/s]（防止差分尖峰灌进指令）
    'commandFilterNaturalFrequency', 20.0, ... % command_filter 的自然频率 [rad/s]
    'commandFilterDampingRatio', 1.0, ...      % command_filter 的阻尼比 [-]
    'highGainDifferentiatorLambda1', 20.0, ...  % 高增益微分器 lambda_1
    'highGainDifferentiatorLambda2', 100.0, ... % 高增益微分器 lambda_2
    'highGainDifferentiatorSmoothing', 1e-3, ... % tanh 平滑尺度 [N]
    'forceDerivativeLimit', 1000.0, ...         % 力导数估计限幅 [N/s]
    'useIntegral', false, ...         % 外层姿态积分开关（与 v2 的 useIntegral 默认一致：关闭）
    'kIR', [0.15; 0.15; 0.10], ...    % 外层姿态积分增益（仅 useIntegral = true 时生效）
    'integralLimit', [0.5; 0.5; 0.5], ... % 外层姿态积分限幅 [rad*s]
    'maxBodyRateCommand', [5.0; 5.0; 3.5], ... % 角速度指令限幅 [rad/s]
    'maxAttitudeError', 2.0, ...      % 姿态误差范数限幅 [rad]
    'headingSource', 'reference');       % 机体航向来源：'reference'(论文原式) | 'worldX'(锁定+x)

% --------------------------- 内环速率跟踪（等效 Crazyflie 内部速率 PID 闭环）
% ★★ 与 v2 同构：**PI**（比例 + 积分），而不是纯比例一阶环节。
%    比例部分保留本项目原值 bandwidth；积分部分取 v2 的 ki 量级。
%    动力学上就是   Omega_dot_i = bandwidth .* e_i + integralGain .* ∫e_i，
%    其中 e_i = Omega_cmd_i - Omega_i（机体惯量不需要，见上面 cfg.vehicle 的说明）。
%    ★ 积分器有限幅；但没有做"抗饱和 gate"。如果观察到指令长期饱和后恢复时
%      出现慢速过冲，参考独立控制器里独立位置积分器抗饱和的做法（加误差门控）。
cfg.rateLoop = struct(...
    'bandwidth', [35.0; 35.0; 20.0], ...  % 比例部分，角速度环带宽 [1/s]
    'integralGain', [3.0; 3.0; 2.0], ...  % 积分部分 [1/s^2]（取 v2 的 ki 量级）
    'integralLimit', [1.5; 1.5; 1.0]);    % 积分限幅 [rad]（v2 的值）

% --------------------------------------------------------------- 初始条件
% 给负载一个明显的初始位置/姿态偏差，各各绳也给不同的初始倾角，
% 用来验证控制器的收敛能力。
cfg.initial = struct();
cfg.initial.position = [0.10; -0.06; ...
    cfg.takeoff.groundZ - 0.5 * cfg.payload.size(3)]; % 负载底面初始接触地面
cfg.initial.velocity = zeros(3, 1);
cfg.initial.rpy = deg2rad([3; -2; 4]);             % 负载初始姿态
cfg.initial.loadBodyRate = zeros(3, 1);

% 初始绳向在 mergeStruct 后根据挂点、绳长和外张间隙自动生成。
% q_i 仍定义为“由四旋翼指向负载”，竖直悬停时 q_i = +e3；允许外张时
% q_i 会带有指向负载中心的水平分量。
cfg.initial.linkUnits = repmat([0; 0; 1], 1, 3);
cfg.initial.linkRates = zeros(3, 3);

% 各机初始姿态与角速度
% ★ vehicleRpy 必须是 3 x n：**行 = [roll; pitch; yaw]，列 = 第 i 架**。
%   曾经写成 deg2rad([5,-4; -3,6; 4,2].')，那是一个 2x3 矩阵，导致第 3 架取
%   出的列向量只有 2 个元素，rpyToRotm 里 rpy(3) 越界。
cfg.initial.vehicleRpy = deg2rad([[5; -4; 4], [-3; 6; 1], [4; 2; -3]]);   % 3 x n，每列一架
cfg.initial.bodyRates = zeros(3, 3);
% 初始推力 = 各机的悬停解。★ 必须按**挂点几何**解，不能假设"三等分"：
%   由 (20)+(17)，悬停时各绳 μ_i 需同时满足力平衡 Σμ = m0 g 与力矩平衡
%   Σ ρ̂_i μ_i = 0。挂点质心若与负载质心重合，三根绳等分；
%   本例新挂点（某边中点 + 对边两顶点）质心偏离 0.0333 m，解出 **2:1:1**：
%       μ = [0.3924, 0.1962, 0.1962] N -> 推力 = |μ| + m g = [0.7112, 0.5150, 0.5150] N
%   若仍用 (m0/3)g + mg = 0.5804 N，机 1 被低估 0.13 N，而机体姿态环立刻要纠偏
%   —— 实测起步瞬间推力打到 **100%**（饱和）。改成按几何解出的向量后降到 70.3%。
rhoInit = cfg.payload.attachPoints;                        % 3 x n
muHoverInit = [ones(1, cfg.vehicle.count); ...
               rhoInit(2, :); ...
               -rhoInit(1, :)] \ ...
              [-cfg.payload.mass * cfg.vehicle.gravity; 0; 0];
cfg.initial.thrustNewton = abs(muHoverInit).' + cfg.vehicle.mass * cfg.vehicle.gravity;

% --------------------------------------------------- 目标点（静态悬停点）
cfg.target = struct();
cfg.target.position = [0; 0; -0.35];               % z 向下为正，-0.35 表示升高
cfg.target.velocity = zeros(3, 1);
cfg.target.acceleration = zeros(3, 1);
cfg.target.rpy = deg2rad([0; 0; 9]);

% ------------------------------------------------------------ 期望参考
% ★ 本工程只有**定高**工况：期望参考就是上面的静态目标点 cfg.target，
%   没有外部轨迹函数（八字轨迹只存在于另一个分支，已连同 cfg.referenceFcn
%   和 cfg.figureEight 一起删除）。
%   两条**时变**参考由仿真主程序自己生成，不在这里配置：
%     ① 交接抬升段：冻结起点 + 5 次多项式剖面（cfg.takeoff.referenceLiftTime）
%     ② 末段下降  ：5 次多项式剖面（cfg.takeoff.landingApproachFraction）
%   两者的 jerk 都是解析算出来的 ⇒ analytic 前馈需要的量依然齐备。
cfg.reference.omegaD = zeros(3, 1);    % 静态目标 ⇒ 负载期望角速度恒为 0
cfg.reference.omegaDotD = zeros(3, 1);

% ------------------------------------------------------------ 可视化参数
%
% ★★★ 负载显示缩放必须为 1.0 —— 这是**几何保真**要求，不是审美偏好。
%   挂点画在**真实**位置（由 payload.size 和 attachFractions 生成），
%   而四旋翼位置由 veh = x0 + R0*rho_i - l*q_i 决定、绳索连到真实挂点。
%   若把负载按 2 倍画，挂点就会落在"板面 1/4 处"，
%   **看起来像挂在板子中间** —— 用户正是据此判断"挂点没在边缘"。
%   曾经这里写 2.0（为了在远景 3D 图里看清外形），已被证实会产生误导。
%   要看大请用 MATLAB 图窗自带的缩放，不要改这个系数。
%   可视化里已加 PayloadScaleNotOne 警告，一旦 ≠ 1 会明确提示该后果。
%
% ★ 注意：本注释块放在 struct(...) **外面**。MATLAB 的 `...` 续行块里，
%   每一行（**包括注释行**）都必须以 `...` 结尾；把注释插进续行块中间
%   会把语句提前截断，报"无效表达式"。这个坑实测踩过。
cfg.visualization = struct(...
    'plot', true, ...
    'animate', true, ...
    'animationStride', 20, ...
    'saveVideo', false, ...
    'videoFile', 'crazyflie_slung_multi_demo.mp4', ...
    'payloadSizeScale', 1.0, ...      % 负载显示尺寸缩放（★必须 1.0 才与挂点一致）
    'vehicleScale', 1.6, ...          % 四旋翼显示尺寸缩放（纯显示，不影响挂点）
    'axisPadding', 0.30, ...          % 三维坐标轴留白 [m]
    'axisSpanMin', 1.80, ...          % 三维坐标轴最小跨度 [m]
    ...                               %   ★ 这个值必须与负载尺寸/绳长同量级。
    ...                               %   曾经是 2.00 m，结果 0.20 m 见方的负载只占画面
    ...                               %   10 %、20 mm 的厚度只占 1 %，负载被压成一片
    ...                               %   看不出外形的薄板。实测 0.55 m 时负载边长占
    ...                               %   36 %、绳索占 64 %，三维图里才真正"看得见"。
    ...                               %   跨度还会被 3.2*负载对角线 与 1.4*绳长抬高
    ...                               %   （见 visualization.m，并对过小的情况发警告）。
    ...                               %   ★ 八字工况下实际跨度由轨迹与障碍物共同决定，
    ...                               %   会远大于 0.55 m，此时该下限自动失效（这没问题：
    ...                               %   跨度大是必然的，负载可读性由负载显示缩放保证）。
    'trajectoryWindow', 25.0, ...     % 三维轨迹子图只显示前 N 秒 [s]
    'trajectoryMinSpan', 0.50);       % 三维轨迹子图的最小轴框跨度 [m]

% ------------------------------------------------------- 参数合法性检查
if cfg.vehicle.count < 3
    error('crazyflie_slung_parameters:NotEnoughVehicles', ...
        ['本仿真实现的是论文的多机协同吊运，需要 n >= 3。' ...
         'n = 2 时 P (6x6) 虽满秩，但存在由 ' ...
         '[(rho_1-rho_2); (rho_2-rho_1)] 张成的一维零空间，' ...
         '张力分配不唯一（会产生内力环流），论文 (2014) 明确规定 n >= 3。' ...
         '当前 cfg.vehicle.count = %d。'], cfg.vehicle.count);
end

% 用户参数覆盖默认值（递归合并）
cfg = mergeStruct(cfg, userCfg);

validOmegaCMethods = {'analytic', 'command_filter', 'high_gain_observer', ...
    'filtered_log_difference', 'none'};
if ~(ischar(cfg.attitudeController.omegaCMethod) ...
        || (isstring(cfg.attitudeController.omegaCMethod) ...
        && isscalar(cfg.attitudeController.omegaCMethod)))
    error('crazyflie_slung_parameters:BadOmegaCMethod', ...
        'attitudeController.omegaCMethod 必须是字符串。');
end
cfg.attitudeController.omegaCMethod = char(cfg.attitudeController.omegaCMethod);
if ~any(strcmpi(cfg.attitudeController.omegaCMethod, validOmegaCMethods))
    error('crazyflie_slung_parameters:BadOmegaCMethod', ...
        '未知 omegaCMethod "%s"；可选值：%s。', ...
        cfg.attitudeController.omegaCMethod, strjoin(validOmegaCMethods, ', '));
end

% ------------------------------------------------------------------ 尺寸驱动的派生参数
% payload.size 是几何参数源。除非用户明确给出更高优先级的自定义值，
% 挂点、惯量、姿态增益和初始推力都在这里重新生成，避免只改尺寸后仍沿用旧系数。
if cfg.vehicle.count < 3
    error('crazyflie_slung_parameters:NotEnoughVehicles', ...
        '本仿真实现的张力分配要求至少 3 架无人机，当前 n = %d。', ...
        cfg.vehicle.count);
end
userHasPayload = isfield(userCfg, 'payload') && isstruct(userCfg.payload);
userHasMass = userHasPayload && isfield(userCfg.payload, 'mass');
userHasDensity = userHasPayload && isfield(userCfg.payload, 'density');
userHasAttachPoints = userHasPayload && isfield(userCfg.payload, 'attachPoints');
userHasInertia = userHasPayload && isfield(userCfg.payload, 'inertia');
if numel(cfg.payload.size) ~= 3 || any(cfg.payload.size(:) <= 0)
    error('crazyflie_slung_parameters:BadPayloadSize', ...
        'cfg.payload.size 必须是三个正数 [length; width; height]。');
end
cfg.payload.size = cfg.payload.size(:);
if userHasDensity && (~isscalar(cfg.payload.density) ...
        || ~isfinite(cfg.payload.density) || cfg.payload.density <= 0)
    error('crazyflie_slung_parameters:BadPayloadDensity', ...
        'cfg.payload.density 必须是正的有限标量 [kg/m^3]。');
end
if userHasDensity && ~userHasMass
    cfg.payload.mass = cfg.payload.density * prod(cfg.payload.size);
end
if ~isscalar(cfg.payload.mass) || ~isfinite(cfg.payload.mass) || cfg.payload.mass <= 0
    error('crazyflie_slung_parameters:BadPayloadMass', ...
        'cfg.payload.mass 必须是正的有限标量 [kg]。');
end
cfg.payload.volume = prod(cfg.payload.size);
cfg.payload.surfaceArea = 2 * (cfg.payload.size(1) * cfg.payload.size(2) ...
    + cfg.payload.size(1) * cfg.payload.size(3) ...
    + cfg.payload.size(2) * cfg.payload.size(3));
cfg.payload.boundingRadius = 0.5 * norm(cfg.payload.size);
if userHasMass || ~userHasDensity || ~isfield(cfg.payload, 'density')
    cfg.payload.density = cfg.payload.mass / cfg.payload.volume;
end
if size(cfg.payload.attachFractions, 1) ~= 3
    error('crazyflie_slung_parameters:BadAttachFractions', ...
        'cfg.payload.attachFractions 必须是 3 x n 矩阵。');
end
if size(cfg.payload.attachFractions, 2) ~= cfg.vehicle.count
    error('crazyflie_slung_parameters:AttachFractionMismatch', ...
        'attachFractions 列数 (%d) 必须等于四旋翼数量 (%d)。', ...
        size(cfg.payload.attachFractions, 2), cfg.vehicle.count);
end
if ~userHasAttachPoints
    cfg.payload.attachPoints = bsxfun(@times, cfg.payload.size(:), ...
        cfg.payload.attachFractions);
end
if ~userHasInertia
    cfg.payload.inertia = boxInertia(cfg.payload.mass, cfg.payload.size);
end
if size(cfg.payload.attachPoints, 1) ~= 3
    error('crazyflie_slung_parameters:BadAttachPoints', ...
        'cfg.payload.attachPoints 必须是 3 x n 矩阵。');
end

% 典型水平力臂 rhoTyp：姿态力矩靠差分张力交付，力臂就是挂点水平投影的最大值。
%   下面两处共用（空中外张的最小倾角地板、姿态增益的力矩预算上限）。
rhoTyp = 0;
for i = 1:size(cfg.payload.attachPoints, 2)
    rhoTyp = max(rhoTyp, norm(cfg.payload.attachPoints(1:2, i)));
end
cfg.payload.horizontalArm = rhoTyp;

% 若用户没有手动指定惯量阻尼，则根据尺寸更新负载的派生物理系数。
if ~isfield(cfg.payload, 'rotationalDampingReferenceSize')
    cfg.payload.rotationalDampingReferenceSize = [0.08; 0.06; 0.050];
end
userHasDamping = userHasPayload && isfield(userCfg.payload, 'rotationalDamping');
if ~userHasDamping
    refSize = cfg.payload.rotationalDampingReferenceSize(:);
    areaScale = (cfg.payload.size(1) * cfg.payload.size(2) ...
        / (refSize(1) * refSize(2)))^2;
    cfg.payload.rotationalDamping = 1.0e-3 * areaScale;
end

% 姿态增益按新的 J0 重算；用户显式给出的 kR/kOmega 保留。
userHasLoadController = isfield(userCfg, 'loadController') ...
    && isstruct(userCfg.loadController);
userHasKR = userHasLoadController && isfield(userCfg.loadController, 'kR');
userHasKOmega = userHasLoadController && isfield(userCfg.loadController, 'kOmega');
wnLoad = 2 * pi * cfg.loadController.designBandwidthHz(:);
zetaLoad = cfg.loadController.designDampingRatio;
manualKR     = cfg.loadController.manualKR;
manualKOmega = cfg.loadController.manualKOmega;
cfg.loadController.attitudeGainManual = ~isempty(manualKR) || ~isempty(manualKOmega) ...
    || userHasKR || userHasKOmega;
if cfg.loadController.attitudeGainManual
    % 手动/用户指定：原样采用，不做自动计算
    if ~isempty(manualKR),     cfg.loadController.kR     = manualKR(:);     end
    if ~isempty(manualKOmega), cfg.loadController.kOmega = manualKOmega(:); end
else
    % 自动：按最终 J0 与目标带宽反解
    cfg.loadController.kR     = diag(cfg.payload.inertia) .* wnLoad.^2;
    cfg.loadController.kOmega = 2 * zetaLoad .* wnLoad .* diag(cfg.payload.inertia);
end

% ======================================================================
% 姿态力矩的『张力预算上限』（2026-09-28 新增）
%   由张力预算反推 kR 上限；超出就削到上限，并按阻尼比重算 kOmega。
%   用户显式给出 kR/kOmega 时不削；小负载不触发
%   （默认尺寸下 kRx = 0.0238 < 上限 0.0567）。
% ======================================================================
momentCap = cfg.loadController.attitudeMomentBudget ...
    * cfg.payload.mass * cfg.vehicle.gravity * rhoTyp;
cfg.loadController.attitudeMomentCap = momentCap;
cfg.loadController.kRCap = momentCap / cfg.loadController.attitudeMomentRefError;
if ~cfg.loadController.attitudeGainManual && any(cfg.loadController.kR > cfg.loadController.kRCap)
    cfg.loadController.kR = min(cfg.loadController.kR, cfg.loadController.kRCap);
    % 保持阻尼比 zeta：kOmega = 2*zeta*sqrt(kR*J0)
    cfg.loadController.kOmega = 2 * zetaLoad ...
        .* sqrt(cfg.loadController.kR .* diag(cfg.payload.inertia));
    cfg.loadController.attitudeGainCapped = true;
else
    cfg.loadController.attitudeGainCapped = false;
end

% 尺寸防火墙：kR/kOmega/kx/kv/ki 必须都是 3x1 列向量。
%   教训（2026-09-28）：曾把 kOmega 写成 sqrt(...).' 得到 1x3，
%   后续 kOmega .* eOmega0 因广播变 3x3 ⇒ 控制器第 125 行
%   rhs6 = [R0.'*Fd; Md] 报『要串联的数组的维度不一致』。
%   该分支只在姿态力矩上限触发时执行，小尺寸不暴露 ⇒ 极难发现。
%   这里显式挡住，让错误在参数阶段就暴露。
gainOK = isequal(size(cfg.loadController.kR), [3, 1]) ...
    && isequal(size(cfg.loadController.kOmega), [3, 1]) ...
    && isequal(size(cfg.loadController.kx), [3, 1]) ...
    && isequal(size(cfg.loadController.kv), [3, 1]) ...
    && isequal(size(cfg.loadController.ki), [3, 1]);
if ~gainOK
    error('crazyflie_slung_parameters:BadGainShape', ...
        ['cfg.loadController 的 kR/kOmega/kx/kv/ki 必须都是 3x1 列向量。' ...
         '当前 kR=%dx%d, kOmega=%dx%d, kx=%dx%d, kv=%dx%d, ki=%dx%d。'], ...
        size(cfg.loadController.kR, 1), size(cfg.loadController.kR, 2), ...
        size(cfg.loadController.kOmega, 1), size(cfg.loadController.kOmega, 2), ...
        size(cfg.loadController.kx, 1), size(cfg.loadController.kx, 2), ...
        size(cfg.loadController.kv, 1), size(cfg.loadController.kv, 2), ...
        size(cfg.loadController.ki, 1), size(cfg.loadController.ki, 2));
end

% 允许无人机在挂点外侧悬挂，而不是强制位于挂点正上方。
if isempty(cfg.link.initialOutwardOffset) || ~isscalar(cfg.link.initialOutwardOffset) ...
        || ~isfinite(cfg.link.initialOutwardOffset)
    cfg.link.initialOutwardOffset = cfg.vehicle.armLength ...
        + cfg.vehicle.rotorRadius + cfg.link.vehicleClearance;
end
% 空中外张 = max(机体包络+间隙, 最小倾角比例 * 最大水平力臂)。
%   原值只与机体有关（尺寸的 0 次方），而力臂正比于尺寸
%   ⇒ 大负载时缆绳相对更竖直、水平张力可用量下降。
%   加最小倾角地板，保证大负载也有同样的倾角比例。
cfg.link.initialOutwardOffset = max(cfg.link.initialOutwardOffset, ...
    cfg.link.minInFlightTiltRatio * rhoTyp);
cfg.link.initialOutwardOffset = min(max(cfg.link.initialOutwardOffset, 0), ...
    0.85 * cfg.link.length);
if ~isfield(cfg.vehicle, 'collisionRadius') || isempty(cfg.vehicle.collisionRadius)
    cfg.vehicle.collisionRadius = cfg.vehicle.armLength + cfg.vehicle.rotorRadius;
end
if ~isfield(cfg.takeoff, 'groundRadialOffset') ...
        || ~isscalar(cfg.takeoff.groundRadialOffset) ...
        || ~isfinite(cfg.takeoff.groundRadialOffset)
    cfg.takeoff.groundRadialOffset = max(cfg.link.initialOutwardOffset, ...
        cfg.payload.boundingRadius + cfg.vehicle.collisionRadius ...
        + cfg.link.vehicleClearance);
end
cfg.takeoff.groundRadialOffset = min(max(cfg.takeoff.groundRadialOffset, 0), ...
    0.85 * cfg.link.length);

% ★★ 收紧段绳向（负载体系，3 x n）—— 由 **simulation.m 在初始化时**填入：
%   它调用控制器自身算一遍悬停平衡解，取那时的期望绳向。
%   这样"收紧段把无人机放到哪"与"控制器期望绳在哪"是**同一个来源**，
%   不会像以前那样由一个统一的 groundRadialOffset 冒充动力学绳向。
%   为空时 simulation.m 会退回旧的统一外张几何（并发 warning）。
cfg.link.takeupLinkUnitsBody = [];

% ======================================================================
% ★★★ 外张机制按『碰撞缺口』自适应 —— 2026-09-27 新增
%
% 动机：`initialOutwardOffset` 与 `outwardBiasFraction` 原本都是**绝对量**，
%   与负载尺寸无关。但『无人机之间够不够开』这件事**只在小负载时才需要外张**：
%     自然机间距 ≈ 挂点间距（随负载尺寸线性增长）
%     要求机间距 = 2*collisionRadius + vehicleClearance（常数）
%   实测本参数集（collisionRadius = 0.0695、clearance = 0.02 ⇒ 要求 0.159 m）：
%     负载边长 0.08 m：无外张时机间距 0.080 m < 0.159 ⇒ 外张**必需**
%     负载边长 0.15 m：0.150 m < 0.159           ⇒ 外张**必需**（临界）
%     负载边长 >= 0.20 m：>= 0.200 m >= 0.159     ⇒ 外张**冗余**
%   而外张内部力的代价是**恒定**的：
%     min(0.20*m0*g/sqrt(n), outwardBiasMax) = 0.0906 N/机
%     = 悬停总张力的 11.5%，但 = **2:1:1 里最小那根绳悬停张力的 46.2%**
%   ⇒ 大负载下它只付代价、不拿收益（白白吃掉最小绳的张力裕度）。
%
% 做法：用『缺口比例』缩放两个外张量；大负载时自动归零。
%   ★ 用户若**显式**给了 outwardBiasFraction / initialOutwardOffset，则尊重用户值，不再缩放。
%   ★ 控制器无需改动：outwardInternalBias 里的 desired 本身按 fraction 缩放，
%     fraction 归零 ⇒ desired 归零 ⇒ bias 归零（其内部 `targetRms > 0` 的守卫不会拦）。
% ======================================================================
reqSep = 2 * cfg.vehicle.collisionRadius + cfg.link.vehicleClearance;
natSep = payloadNaturalSeparation(cfg.payload.attachPoints);
cfg.allocation.collisionRequiredSeparation = reqSep;
cfg.allocation.payloadNaturalSeparation  = natSep;
cfg.allocation.outwardBiasScale = min(max((reqSep - natSep) / max(reqSep, eps), 0), 1);

userHasAlloc = isfield(userCfg, 'allocation') && isstruct(userCfg.allocation);
if ~(userHasAlloc && isfield(userCfg.allocation, 'outwardBiasFraction'))
    cfg.allocation.outwardBiasFraction = cfg.allocation.outwardBiasFraction ...
        * cfg.allocation.outwardBiasScale;
end
if ~(isfield(userCfg, 'link') && isstruct(userCfg.link) ...
        && isfield(userCfg.link, 'initialOutwardOffset'))
    cfg.link.initialOutwardOffset = cfg.link.initialOutwardOffset ...
        * cfg.allocation.outwardBiasScale;
end
userHasInitial = isfield(userCfg, 'initial') && isstruct(userCfg.initial);
userHasInitialPosition = userHasInitial && isfield(userCfg.initial, 'position');
userHasInitialLinks = userHasInitial && isfield(userCfg.initial, 'linkUnits');
userHasInitialRates = userHasInitial && isfield(userCfg.initial, 'linkRates');
userHasInitialBodyRates = userHasInitial && isfield(userCfg.initial, 'bodyRates');
userHasInitialThrust = userHasInitial && isfield(userCfg.initial, 'thrustNewton');

% 起飞模型的派生几何和数值检查。若用户没有显式给 initial.position，
% 始终把负载底面放在 groundZ，避免修改负载厚度后初始位置悬空或穿地。
if ~userHasInitialPosition
    cfg.initial.position = [cfg.initial.position(1); cfg.initial.position(2); ...
        cfg.takeoff.groundZ - 0.5 * cfg.payload.size(3)];
end
if ~isfield(cfg, 'takeoff') || ~isstruct(cfg.takeoff)
    error('crazyflie_slung_parameters:BadTakeoffConfig', ...
        'cfg.takeoff 必须是结构体。');
end
if cfg.takeoff.epsilonOff <= cfg.takeoff.epsilonOn
    error('crazyflie_slung_parameters:BadTakeoffHysteresis', ...
        'takeoff.epsilonOff 必须大于 epsilonOn。');
end
% ★★ preTensionSlack 必须允许 0（默认就是 0），理由见下面的说明 ★★
%   本模型**只能表示绷紧的绳**（松弛段未实现）。所以"收紧末端保留 σ 的绳长余量"
%   在数学上等价于"进入绷紧段时把 σ 一次性收掉"—— 实测那是一次 **17 mm 的
%   无人机位置瞬移**（等价速率 8.5 m/s），直接激励绳向环。
%   ⇒ 正确的交接点是"绳刚好拉直、张力为零"，即 σ = 0。
%   （σ > 0 仍然允许，但会在交接瞬间产生同样大小的瞬移，会发 warning。）
if cfg.takeoff.preTensionSlack < 0 || ...
        cfg.takeoff.preTensionSlack >= cfg.link.length
    error('crazyflie_slung_parameters:BadTakeupSlack', ...
        'takeoff.preTensionSlack 必须在 [0, link.length) 内（本模型建议取 0）。');
end

% 用户只改 rpy 时自动重算旋转矩阵
if ~isfield(userCfg, 'initial') || ~isfield(userCfg.initial, 'rpy')
    cfg.initial.R0 = rpyToRotm(cfg.initial.rpy);
elseif ~isfield(userCfg.initial, 'R0')
    cfg.initial.R0 = rpyToRotm(cfg.initial.rpy);
end
if ~isfield(cfg.initial, 'R0')
    cfg.initial.R0 = rpyToRotm(cfg.initial.rpy);
end
% 初始姿态可能使长方体的竖直包络大于 H/2。未显式指定初始位置时，
% 按当前姿态的竖直包络把负载底部放在地面上，避免角部穿地。
if cfg.takeoff.enabled && ~userHasInitialPosition
    halfHeight = 0.5 * sum(abs(cfg.initial.R0(3, :)) ...
        .* cfg.payload.size(:).');
    cfg.initial.position(3) = cfg.takeoff.groundZ - halfHeight;
end
if ~userHasInitialLinks
    if cfg.link.allowTiltedCables
        cfg.initial.linkUnits = defaultLinkUnits(cfg.payload.attachPoints, ...
            cfg.link.length, cfg.link.initialOutwardOffset, cfg.initial.R0);
    else
        cfg.initial.linkUnits = repmat([0; 0; 1], 1, cfg.vehicle.count);
    end
end
if ~userHasInitialRates
    cfg.initial.linkRates = zeros(3, cfg.vehicle.count);
end
if ~userHasInitialBodyRates
    cfg.initial.bodyRates = zeros(3, cfg.vehicle.count);
end
cfg.target.R0 = rpyToRotm(cfg.target.rpy);

% 各机初始旋转矩阵（3x3xn）
if ~isfield(cfg.initial, 'vehicleR') || isempty(cfg.initial.vehicleR)
    if ~isequal(size(cfg.initial.vehicleRpy), [3, cfg.vehicle.count])
        error('crazyflie_slung_parameters:BadVehicleRpy', ...
            ['cfg.initial.vehicleRpy 必须是 3 x n 矩阵（行 = [roll;pitch;yaw]，' ...
             '列 = 第 i 架），当前尺寸为 %dx%d，而 n = %d。'], ...
            size(cfg.initial.vehicleRpy, 1), size(cfg.initial.vehicleRpy, 2), ...
            cfg.vehicle.count);
    end
    cfg.initial.vehicleR = zeros(3, 3, cfg.vehicle.count);
    for i = 1:cfg.vehicle.count
        cfg.initial.vehicleR(:, :, i) = rpyToRotm(cfg.initial.vehicleRpy(:, i));
    end
end

% 挂点数必须与车辆数一致
if size(cfg.payload.attachPoints, 2) ~= cfg.vehicle.count
    error('crazyflie_slung_parameters:AttachPointMismatch', ...
        '挂点数 (%d) 必须等于四旋翼数量 (%d)。', ...
        size(cfg.payload.attachPoints, 2), cfg.vehicle.count);
end
if ~userHasInitialThrust
    rhoInit = cfg.payload.attachPoints;
    muHoverInit = [ones(1, cfg.vehicle.count); ...
                   rhoInit(2, :); ...
                  -rhoInit(1, :)] \ ...
                   [-cfg.payload.mass * cfg.vehicle.gravity; 0; 0];
    % 初始推力按当前 q_i 的向量方向估计。q_i = +e3 时退化为
    % |mu_i| + m_i*g；允许外张时会自动包含所需的水平推力分量。
    e3 = [0; 0; 1];
    qInit = cfg.initial.linkUnits;
    thrustInit = zeros(1, cfg.vehicle.count);
    for i = 1:cfg.vehicle.count
        qi = normalizeVector(qInit(:, i));
        thrustVector = -abs(muHoverInit(i)) * qi ...
            - cfg.vehicle.mass * cfg.vehicle.gravity * e3;
        thrustInit(i) = norm(thrustVector);
    end
    cfg.initial.thrustNewton = thrustInit;
end
cfg.link.count = cfg.vehicle.count;
end

% ======================================================================
% 局部工具函数（详见 README §1.1：MATLAB 局部函数是文件私有的）
% ======================================================================
function J = boxInertia(mass, sizeXYZ)
%BOXINERTIA 均质长方体绕质心的惯量矩阵。
sizeXYZ = sizeXYZ(:);
a = sizeXYZ(1); b = sizeXYZ(2); c = sizeXYZ(3);
J = diag([mass * (b^2 + c^2) / 12; ...
          mass * (a^2 + c^2) / 12; ...
          mass * (a^2 + b^2) / 12]);
end

function value = payloadNaturalSeparation(rhoAll)
%PAYLOADNATURALSEPARATION 三机在『挂点正上方』(不外加外张偏移) 时的最小两两间距。
% 用于判断碰撞约束是否真的需要外张：若该值已 >= 要求间距，则外张纯属多余。
n = size(rhoAll, 2);
value = inf;
for i = 1:n
    for j = i + 1:n
        value = min(value, norm(rhoAll(:, i) - rhoAll(:, j)));
    end
end
if ~isfinite(value)
    value = 0;
end
end

% ======================================================================
function qAll = defaultLinkUnits(rhoAll, linkLength, outwardOffset, R0)
%DEFAULTLINKUNITS 依据挂点径向方向生成外张的初始绳向。
% q_i 从无人机指向负载；无人机向挂点外侧偏移时，q_i 的水平分量指向内侧。
n = size(rhoAll, 2);
qAll = zeros(3, n);
for i = 1:n
    radialBody = [rhoAll(1, i); rhoAll(2, i); 0];
    radialNorm = norm(radialBody);
    if radialNorm < 1e-12
        angle = 2 * pi * (i - 1) / max(n, 1);
        radialBody = [cos(angle); sin(angle); 0];
    else
        radialBody = radialBody / radialNorm;
    end
    horizontalRatio = min(outwardOffset / max(linkLength, eps), 0.85);
    qBody = [-horizontalRatio * radialBody(1); ...
             -horizontalRatio * radialBody(2); ...
              sqrt(max(1 - horizontalRatio^2, 0))];
    qAll(:, i) = normalizeVector(R0 * qBody);
end
end

function R = rpyToRotm(rpy)
% ZYX 欧拉角顺序，与 v2 保持一致的接口。
roll = rpy(1); pitch = rpy(2); yaw = rpy(3);
cr = cos(roll);  sr = sin(roll);
cp = cos(pitch); sp = sin(pitch);
cy = cos(yaw);   sy = sin(yaw);
R = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr; ...
     sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr; ...
     -sp,   cp*sr,            cp*cr];
end

function v = normalizeVector(v)
v = v(:);
n = norm(v);
if n < eps
    v = [0; 0; 1];      % 退化时取 +e3：四旋翼在负载上方（悬停构型）
else
    v = v / n;
end
end

function out = mergeStruct(base, override)
% 递归合并结构体，便于只覆盖少数参数。
out = base;
fields = fieldnames(override);
for i = 1:numel(fields)
    name = fields{i};
    if isstruct(override.(name)) && isfield(base, name) && isstruct(base.(name))
        out.(name) = mergeStruct(base.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end
end
