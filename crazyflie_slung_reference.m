function [position, velocity, acceleration, rotation, bodyRate, bodyRateDot] = ...
    crazyflie_slung_reference(t, cfg)
%CRAZYFLIE_SLUNG_REFERENCE 三阶段八字避障轨迹参考（起飞 -> 环绕两个八字 -> 降落）。
%
% 对应论文 Lee 2018 第 V 节的数值算例（Fig. 3）：
%   "three quadrotors (n = 3) transport a rectangular box along a figure-eight
%    curve around two obstacles represented by cones"
%   期望轨迹    x0d(t) = [1.2 sin(0.2 pi t), 4.2 cos(0.1 pi t), -0.5]'
%   期望姿态    R0d(t) = [ x0d_dot/||x0d_dot|| , hat(e3) x0d_dot/||hat(e3) x0d_dot|| , e3 ]
%               "its first axis is tangent to the desired path, and the third
%                axis is parallel to the direction of gravity"
%
% ★ 本文件相对论文做了四处工程化改造，参数都在 cfg.figureEight 里可调：
%
%   1. **尺度缩放**：论文的八字横宽 ±1.2 m、纵向 ±4.2 m，配负载 1.0x0.8x0.2 m、
%      绳长 1 m。本仿真负载只有 0.20 m、绳长 0.35 m，照搬论文尺度会让负载/绳索
%      在画面里缩成一个点，摆动也会被极端激发。故按 cfg.figureEight 的
%      amplitudeX / amplitudeY 取横宽 ±1.5 m 量级（用户选定"大尺度"）。
%
%   2. **三阶段**（用户要求：起飞 + 环绕两个八字 + 降落）：
%         阶段 0  起飞    t ∈ [0, T1)              竖直爬升到环绕平面
%         阶段 1  环绕    t ∈ [T1, T1+T2)           绕两个完整八字
%         阶段 2  降落    t ∈ [T1+T2, T1+T2+T3]     水平归位 + 竖直下降到落点
%      三段在拼接点上**位置/速度/加速度全部连续**（C^2），见下面的
%      smoothstep 与 blendEnvelope 说明。
%
%   3. **期望偏航角恒定**：论文的 R0d 让负载第一轴指向速度方向，即偏航随路径转。
%      但本仿真已确认负载偏航通道**结构性不可控**（详见 parameters.m 与
%      README §5.1）。若把期望偏航设成时变量，偏航误差会持续激励那个不可控
%      通道并导致发散。故这里把 R0d 的第一轴固定为 +x（期望偏航恒为 0），
%      与静态工况的处理一致。cfg.figureEight.lockYaw = false 可恢复论文原式
%      （仅供研究，实测会发散）。
%
%   4. **环绕段起止点归零**：论文的八字是无限绕行的，但本工况要求"绕两个八字
%      之后降落"。为了让阶段 1 与阶段 2 的衔接点在**同一水平位置**，
%      这里把八字参数化平移使 τ = 0 与 τ = tEnd 都落在原点：
%          pX_bar(τ) = aX * sin(wX τ)                  (τ=0 与 τ=tEnd 均为 0)
%          pY_bar(τ) = aY * (1 - cos(wY τ)) / 2        (τ=0 与 τ=tEnd 均为 0)
%      注意这是**保形变换**：曲线仍然是标准双纽线（Lemniscate of Gerono），
%      只是把坐标系原点移到了自交点、并把 y 幅值从 aY 缩到 aY/2。
%      ★ 论文原式 x = aX sin(wX t), y = aY cos(wY t) 的自交点不在 t=0，
%        所以"绕完两个八字回到出发点"这个要求必须做这个平移。
%
% 输出（与 crazyflie_slung_simulation.m 的 referenceState 约定一致）：
%   position     (3x1)  x0d
%   velocity     (3x1)  x0d_dot
%   acceleration (3x1)  x0d_ddot
%   rotation     (3x3)  R0d
%   bodyRate     (3x1)  Omega0d（负载体系角速度）
%   bodyRateDot  (3x1)  Omega0d_dot（负载体系角加速度）
%
% 用法（parameters.m 里已自动配置）：
%   cfg.referenceFcn = @(t) crazyflie_slung_reference(t, cfg);
%   sim = crazyflie_slung_simulation(cfg);

% ------------------------------------------------------------------ 参数
fe = cfg.figureEight;
T1 = fe.takeoffDuration;      % 起飞段时长 [s]
T2 = fe.cruiseDuration;       % 环绕段时长 [s]（覆盖 fe.cycles 个完整八字）
T3 = fe.landingDuration;      % 降落段时长 [s]

aX = fe.amplitudeX;           % 半宽（x 方向）[m]
aY = fe.amplitudeY;           % 半长（y 方向）[m]
wX = fe.omegaX;               % x 角频率 [rad/s]，论文为 0.2*pi
wY = fe.omegaY;               % y 角频率 [rad/s]，论文为 0.1*pi

zCruise = fe.cruiseHeight;    % 环绕平面高度（z 向下为正，故为负值）
xyStart = fe.startPosition(1:2);   % 起飞起点水平位置
zStart = fe.startPosition(3);
xyLand = fe.landPosition(1:2);     % 降落终点水平位置
zLand = fe.landPosition(3);

% 八字段的定义域：让环绕段正好覆盖 fe.cycles 个完整周期。
% y 的一个完整周期是 2*pi/wY（论文为 20 s），故 tEnd = cycles * 2*pi/wY。
tEnd = fe.cycles * 2 * pi / wY;

% ------------------------------------------------------------------ 相位
if t < T1
    phase = 0;                  % 起飞
elseif t < T1 + T2
    phase = 1;                  % 环绕
else
    phase = 2;                  % 降落
end

% ------------------------------------------------ 八字轨迹（平移+缩放后的保形形式）
% 时间映射：环绕段把 [T1, T1+T2] 线性映射到 [0, tEnd]
%   τ = (t - T1) * s,  s = tEnd / T2
% 因此对**真实时间**求导要乘 s（速度）、s^2（加速度）——链式法则。
s = tEnd / T2;
if phase == 0
    tau = 0;
elseif phase == 1
    tau = (t - T1) * s;
else
    tau = tEnd;
end

% 平移后的八字（τ = 0 与 τ = tEnd 都经过原点）
%   pX_bar =  aX * sin(wX τ)              pX_bar(0) = 0, pX_bar(tEnd) = 0
%   pY_bar =  aY * (1 - cos(wY τ)) / 2    pY_bar(0) = 0, pY_bar(tEnd) = 0
pHatX = aX * sin(wX * tau);
pHatY = 0.5 * aY * (1 - cos(wY * tau));
vHatX = aX * wX * cos(wX * tau);
vHatY = 0.5 * aY * wY * sin(wY * tau);
cHatX = -aX * wX^2 * sin(wX * tau);
cHatY = 0.5 * aY * wY^2 * cos(wY * tau);

% ------------------------------------------------------------------ 三阶段
switch phase
    case 0
        % ---------------- 阶段 0：起飞 ----------------
        % 从 startPosition 平滑走到**环绕段起点**（即平移后八字在 τ=0 的点，
        % 也就是 (0, 0, zCruise)）。用 5 次 smoothstep 保证两端
        % 位置/速度/加速度连续，且首尾速度与加速度都为 0。
        %   sigma(u)  = 6u^5 - 15u^4 + 10u^3   -> sigma(0)=0, sigma(1)=1
        %   sigma'(u) = 30u^2 (1-u)^2          -> sigma'(0)=sigma'(1)=0
        %   sigma''(u)= 60u(1-u)(1-2u)         -> sigma''(0)=sigma''(1)=0
        u = t / T1;
        [sigma, sigmaDot, sigmaDDot] = smoothstep(u);
        targetXY = [0; 0];                 % 环绕段起点（平移后八字在 τ=0）
        position = [xyStart(1) + (targetXY(1) - xyStart(1)) * sigma; ...
                    xyStart(2) + (targetXY(2) - xyStart(2)) * sigma; ...
                    zStart + (zCruise - zStart) * sigma];
        velocity = [(targetXY(1) - xyStart(1)) * sigmaDot / T1; ...
                    (targetXY(2) - xyStart(2)) * sigmaDot / T1; ...
                    (zCruise - zStart) * sigmaDot / T1];
        acceleration = [(targetXY(1) - xyStart(1)) * sigmaDDot / T1^2; ...
                        (targetXY(2) - xyStart(2)) * sigmaDDot / T1^2; ...
                        (zCruise - zStart) * sigmaDDot / T1^2];

    case 1
        % ---------------- 阶段 1：环绕两个八字 ----------------
        % ★ 为什么需要 blendEnvelope：
        %   八字是**闭环**曲线，在 τ=tEnd 处的速度 vHatX(tEnd) = aX*wX ≠ 0；
        %   而起飞段末速度必须为 0（smoothstep 性质）、降落段初速度也必须是 0。
        %   若直接把八字接上去，两个拼接点都会出现速度跳变（实测 2.09 m/s），
        %   控制器会看到阶跃加速度指令而剧烈抖振。
        %   故用**两端对称**的速度包络 env(τ_local)：
        %     前 blendTime 秒由 0 升到 1，末尾 blendTime 秒由 1 降到 0，
        %     且 env、env'、env'' 在各段端点都为 0。
        %   于是两个拼接点处
        %     位置 = env*p = 0            与相邻段末/初位置一致（都是原点）
        %     速度 = env'*p + env*p'*s = 0 与相邻段末/初速度一致
        %     加速度同理为 0               与相邻段末/初加速度一致
        %   blendTime 之外的中段 env ≡ 1，轨迹与论文的八字完全一致。
        [env, envDot, envDDot] = blendEnvelope(t - T1, T2, fe.blendTime);
        position = [pHatX * env; pHatY * env; zCruise];
        % pos = env * p(τ),  τ = (t-T1)*s  =>  d/dt = env' * p + env * p' * s
        %   d2/dt2 = env'' * p + 2 env' p' s + env p'' s^2   （乘积法则）
        velocity = [envDot * pHatX + env * s * vHatX; ...
                    envDot * pHatY + env * s * vHatY; ...
                    0];
        acceleration = [envDDot * pHatX + 2 * envDot * s * vHatX + env * s^2 * cHatX; ...
                        envDDot * pHatY + 2 * envDot * s * vHatY + env * s^2 * cHatY; ...
                        0];

    otherwise
        % ---------------- 阶段 2：降落 ----------------
        % 从环绕段终点（τ=tEnd，平移后八字在 (0,0,zCruise)）平滑走到
        % landPosition。水平与竖直可以各自用 smoothstep 一起插值 ——
        % 因为起点与终点都在同一平面高度上，水平归位与降高同时进行，
        % 视觉上像"飞回落点并降落"。
        % 两端的速度与加速度都为 0，与环绕段末态（env≡1、但八字在 tEnd
        % 处速度 vHatX(tEnd) = aX wX cos(wX tEnd) ≠ 0）—— 见下方说明。
        u = (t - T1 - T2) / T3;
        [sigma, sigmaDot, sigmaDDot] = smoothstep(u);
        fromXY = [0; 0];
        position = [fromXY(1) + (xyLand(1) - fromXY(1)) * sigma; ...
                    fromXY(2) + (xyLand(2) - fromXY(2)) * sigma; ...
                    zCruise + (zLand - zCruise) * sigma];
        velocity = [(xyLand(1) - fromXY(1)) * sigmaDot / T3; ...
                    (xyLand(2) - fromXY(2)) * sigmaDot / T3; ...
                    (zLand - zCruise) * sigmaDot / T3];
        acceleration = [(xyLand(1) - fromXY(1)) * sigmaDDot / T3^2; ...
                        (xyLand(2) - fromXY(2)) * sigmaDDot / T3^2; ...
                        (zLand - zCruise) * sigmaDDot / T3^2];
end

% ------------------------------------------------ 期望姿态（论文 R0d 构造）
% 论文原式：第一轴 = **运动方向**，第三轴 = 重力方向 e3。
%
% ★★★ 巡航段期望偏航 = **参考轨迹的实际运动方向**，|v|→0 处平滑过渡。
%   为什么不用"理想八字切向"：实际 v = env'·p_hat + env·s·v_hat，在包络过渡段
%   |v|→0 时被 **env'·p_hat** 主导 ⇒ 真实运动是**径向**、理想切向是**切向**
%   ⇒ 实测最大差 93.1° / 162.7° / 178.3°（τ=3.0/43.0/44.5 s），
%     表现就是"偏航有时朝运动方向、有时正好相反"。
%   为什么不能直接 atan2(vy,vx)：|v|→0 处方向会翻转（实测 t=47.982 s 处
%   psi 由 +179.98° 一步跳到 0°，|v| 仅 1e-6）。
%   做法：w = v + epsV*t_hat_unit（t_hat = 夹紧 tau 后的理想切向，恒不为零），
%         psi = atan2(w_y, w_x)：|v|>>epsV ⇒ 就是实际运动方向；
%                                 |v|→0    ⇒ 平滑接过渡方向，不跳 180°。
%   实测 |v|>0.02 时与真实运动方向最大差 **0.147°**。
%
% ★ 起飞/降落段仍取 psi = 0（不跟随）：① 近零速、方向病态；
%   ② 巡航两端 |v|→0 的极限方向**恰为 +x**，与 psi=0 天然衔接 ⇒ 无跳变。
%   ⚠ 曾推广到全阶段 ⇒ 起飞段水平方向 149° 让参考 yaw 转 149° ⇒ 拧绳 ⇒
%     最小张力裕度 33.0%→1.9%、推力峰值 69.1%→100%。⇒ 只保留巡航段。
if fe.lockYaw || phase ~= 1
    psiYaw = 0;  psiDot = 0;  psiDotDot = 0;
else
    tauClamped = min(max(tau, 0), fe.cycles * 2 * pi / fe.omegaY);
    tHatX = aX * wX * cos(wX * tauClamped);
    tHatY = 0.5 * aY * wY * sin(wY * tauClamped);
    tNorm = hypot(tHatX, tHatY);
    if tNorm > 0
        tHatX = tHatX / tNorm;  tHatY = tHatY / tNorm;
    end
    epsV = 1e-4;
    wx = velocity(1) + epsV * tHatX;
    wy = velocity(2) + epsV * tHatY;
    den = wx^2 + wy^2;
    psiYaw = atan2(wy, wx);
    if den > 1e-18
        psiDot = (wx * acceleration(2) - wy * acceleration(1)) / den;
    else
        psiDot = 0;
    end
    h = max(cfg.simulation.dt, 1e-4);
    psiDotDot = (linkYawRateAt(t + h, T1, T2, fe) ...
               - linkYawRateAt(t - h, T1, T2, fe)) / (2 * h);
end
b1d = [cos(psiYaw); sin(psiYaw); 0];
b2d = [-sin(psiYaw); cos(psiYaw); 0];
b3d = [0; 0; 1];
rotation = [b1d, b2d, b3d];
bodyRate = [0; 0; psiDot];
bodyRateDot = [0; 0; psiDotDot];

position = position(:);
velocity = velocity(:);
acceleration = acceleration(:);
bodyRate = bodyRate(:);
bodyRateDot = bodyRateDot(:);
end

% ======================================================================
function [vxy, axy] = cruiseHorizontalState(t, T1, T2, fe)
%CRUISEHORIZONTALSTATE 环绕段在真实时间 t 处的**水平**速度与加速度（惯性系）。
% ★ 供 linkYawRateAt 求 psiDotDot 的中心差分；与主函数 case 1 **必须逐式一致**。
aX = fe.amplitudeX;  aY = fe.amplitudeY;
wX = fe.omegaX;      wY = fe.omegaY;
s = fe.cycles * 2 * pi / wY / T2;
tau = (t - T1) * s;
pHatX = aX * sin(wX * tau);           pHatY = 0.5 * aY * (1 - cos(wY * tau));
vHatX = aX * wX * cos(wX * tau);      vHatY = 0.5 * aY * wY * sin(wY * tau);
cHatX = -aX * wX^2 * sin(wX * tau);   cHatY = 0.5 * aY * wY^2 * cos(wY * tau);
[env, envDot, envDDot] = blendEnvelope(t - T1, T2, fe.blendTime);
vxy = [envDot * pHatX + env * s * vHatX; envDot * pHatY + env * s * vHatY];
axy = [envDDot * pHatX + 2 * envDot * s * vHatX + env * s^2 * cHatX; ...
       envDDot * pHatY + 2 * envDot * s * vHatY + env * s^2 * cHatY];
end

% ======================================================================
function d = linkYawRateAt(t, T1, T2, fe)
%LINKYAWRATEAT 环绕段 psiDot（解析式），供其中心差分求 psiDotDot。
[v, a] = cruiseHorizontalState(t, T1, T2, fe);
den = v(1)^2 + v(2)^2;
if den > 1e-12
    d = (v(1) * a(2) - v(2) * a(1)) / den;
else
    d = 0;
end
end

% ======================================================================
function [psi, psiDot, psiDotDot] = yawFromIdealTangent(tau, s, h, aX, aY, wX, wY)
% 理想八字（不含包络）的切向偏航角及其一/二阶**真实时间**导数。
%
%   psi       = atan2(v, u)，  u = aX wX cos(wX tau)， v = (aY/2) wY sin(wY tau)
%   dpsi/dtau = (u v' - v u') / (u^2 + v^2)      （u'、v' 对 tau 求导，解析）
%   psiDot    = dpsi/dtau * s                    （链式法则 d/dt = s d/dtau）
%   psiDotDot = 对**解析的** psiDot 再做一次中心差分
%               （只一阶；比"对 R0d 做二阶差分"的噪声小一个量级）
[psi, psiDot] = yawPsiAndRate(tau, s, aX, aY, wX, wY);
[~, pPlus]  = yawPsiAndRate(tau + s * h, s, aX, aY, wX, wY);
[~, pMinus] = yawPsiAndRate(tau - s * h, s, aX, aY, wX, wY);
psiDotDot = (pPlus - pMinus) / (2 * h);
end

% ======================================================================
function [psi, psiDot] = yawPsiAndRate(tau, s, aX, aY, wX, wY)
% 上面那个函数的解析部分（中心差分需要重复调用它）。
u  = aX * wX * cos(wX * tau);
v  = 0.5 * aY * wY * sin(wY * tau);
uP = -aX * wX^2 * sin(wX * tau);
vP = 0.5 * aY * wY^2 * cos(wY * tau);
psi = atan2(v, u);
psiDot = ((u * vP - v * uP) / max(u^2 + v^2, eps)) * s;
end

% ======================================================================
function [sigma, sigmaDot, sigmaDDot] = smoothstep(u)
% 5 次平滑阶跃：两端的位置、速度、加速度都连续。
% 对 u 超出 [0,1] 的情况做钳位，避免末端数值越界。
u = min(max(u, 0), 1);
sigma = u^3 * (10 - 15 * u + 6 * u^2);
sigmaDot = 30 * u^2 * (1 - u)^2;
sigmaDDot = 60 * u * (1 - u) * (1 - 2 * u);
end

% ======================================================================
function [env, envDot, envDDot] = blendEnvelope(tauLocal, cruiseDuration, blendTime)
% 环绕段的**两端对称**速度包络：0 -> 1 -> 0。
%
% 用途：八字是闭环曲线，在定义域两端速度都不为 0，而相邻的起飞/降落段
% 两端速度必须为 0，否则拼接点出现速度跳变（实测 2.09 m/s）。
% 本包络让八字速度在前 blendTime 秒由 0 升到 1、在末尾 blendTime 秒
% 由 1 降到 0，且 env、env'、env'' 在各段端点都为 0，从而保证
% 两个拼接点的位置、速度、加速度全部连续（C^2）。
%
%   tauLocal       相对环绕段起点的时间 [s]
%   cruiseDuration 环绕段总时长 [s]
%   blendTime      单侧过渡时长 [s]
%
% 形状：env = sigma(u_in) * sigma(u_out)
%   u_in  = tauLocal / blendTime
%   u_out = (cruiseDuration - tauLocal) / blendTime
% sigma 是 5 次 smoothstep。乘积在 u_in -> 0 或 u_out -> 0 时都趋于 0，
% 且两端一阶二阶导都为 0，满足 C^2 拼接要求。
if blendTime <= 0
    env = 1;
    envDot = 0;
    envDDot = 0;
    return;
end
uIn = tauLocal / blendTime;
uOut = (cruiseDuration - tauLocal) / blendTime;
[sIn, sInDot, sInDDot] = smoothstep(uIn);
[sOut, sOutDot, sOutDDot] = smoothstep(uOut);
% 链式法则：d/dtauLocal = (1/blendTime) d/du
env = sIn * sOut;
envDot = (sInDot * sOut + sIn * sOutDot) / blendTime;
envDDot = (sInDDot * sOut + 2 * sInDot * sOutDot + sIn * sOutDDot) / blendTime^2;
end

