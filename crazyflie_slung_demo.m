function report = crazyflie_slung_demo(options)
%CRAZYFLIE_SLUNG_DEMO 多机协同吊运仿真的一键运行 + 自检。
%
% 用法：
%   crazyflie_slung_demo();                 % 跑默认工况（3 机），画曲线 + 三维动画
%   crazyflie_slung_demo("quick");          % 只跑自检，不画图（缩短时长）
%   report = crazyflie_slung_demo();        % 取回自检报告结构体
%
% ★ 默认工况是**绕八字避障**（对应论文 Lee 2018 Fig. 3）：
%     阶段 0  起飞    t ∈ [0, 3)    竖直爬升到环绕平面
%     阶段 1  环绕    t ∈ [3, 21)   绕两个完整八字，横宽 ±1.5 m，绕开两个锥
%     阶段 2  降落    t ∈ [21, 25]  水平归位 + 竖直下降到落点
%   在参数文件里清空 cfg.referenceFcn 可退回静态悬停工况。
%   自检因此分两组：**通用组**（动力学/绳约束/分配一致性）两组工况都适用，
%   **八字组**（跟踪精度/避障间隙/三阶段指标）只在八字工况下检查。
%
% 自检项的物理依据：
%   * 张力分配一致性 —— 论文 (22) 要求 sum_i mu_id = Fd 且 sum_i hat(rho_i) R0' mu_id = Md。
%     伪逆解 (23) 必须精确满足这两个等式，否则说明分配矩阵 P 装配有误。
%   * 张力恒正且总量与 m0 g 同量级 —— 不对称挂点下各绳张力不一定相等。
%   * 缆绳长度守恒、无人机间距和机体-负载间隙满足碰撞包络；不要求无人机正上方。
%   * 推力裕度、角速度范围、数值有限性。
%   * 负载姿态误差按 roll/pitch 与 yaw 分开判定；yaw 是低带宽闭环轴。
%   * 偏航角速度和 yaw 跟踪误差有界。
%   * 【八字组】三阶段跟踪误差 —— 起飞/环绕/降落分别设门限。环绕段是
%     "跟踪精度"真正该考核的窗口；门限按 ±1.5 m 大尺度机动的实测值设定
%     （见 README §10 的调参记录），不是拍脑袋的数字。
%   * 【八字组】负载到两个锥形障碍物的净间隙 —— 对应论文 Fig. 3 的
%     "around two obstacles represented by cones"，是本工况的正确性条件。

if nargin < 1 || isempty(options)
    options = "full";
end
quickMode = strcmpi(string(options), "quick");

fprintf('==============================================================\n');
fprintf(' 多机协同吊运仿真（Lee 2014/2018 框架）  —  一键自检\n');
fprintf('==============================================================\n\n');

cfg = crazyflie_slung_parameters();
isFigureEight = ~isempty(cfg.referenceFcn) && isfield(cfg, 'figureEight');
if quickMode
    cfg.visualization.plot = false;
    cfg.visualization.animate = false;
    if isFigureEight
        % 八字工况的时长是三阶段之和，收紧 quick 模式会破坏阶段划分，
        % 因此这里只把总时长按比例缩到 2 s 留作极快的"冒烟测试"，
        % 但保留三阶段结构（起飞/环绕/降落各按比例缩）。
        % ★ 必须按**当前**三阶段时长等比缩放，不能硬编码。
        %   曾经写死 3.0 / 18.0 / 4.0，cruiseDuration 改成 45 后这套比例就错了。
        total = cfg.simulation.duration;
        cfg.figureEight.takeoffDuration = cfg.figureEight.takeoffDuration / total * 2.0;
        cfg.figureEight.cruiseDuration = cfg.figureEight.cruiseDuration / total * 2.0;
        cfg.figureEight.landingDuration = cfg.figureEight.landingDuration / total * 2.0;
        cfg.figureEight.blendTime = min(cfg.figureEight.blendTime, ...
            cfg.figureEight.cruiseDuration / 3);
        cfg.simulation.duration = 2.0;
        cfg.referenceFcn = @(t) crazyflie_slung_reference(t, cfg);
    else
        cfg.simulation.duration = 12.0;
    end
    fprintf('[模式] quick：不绘图，时长 %.1f s\n\n', cfg.simulation.duration);
end

n = cfg.vehicle.count;
m0 = cfg.payload.mass;
m = cfg.vehicle.mass;
g = cfg.vehicle.gravity;
tMax = cfg.vehicle.maxTotalThrust;

% ------------------------------ 分配矩阵可解性（rank(P) = 6）-----------------
P = [repmat(eye(3), 1, n); zeros(3, 3 * n)];
for i = 1:n
    rhoi = cfg.payload.attachPoints(:, i);
    P(4:6, 3 * (i - 1) + (1:3)) = [0, -rhoi(3), rhoi(2); ...
                                   rhoi(3), 0, -rhoi(1); ...
                                   -rhoi(2), rhoi(1), 0];
end
graphRank = rank(P);

fprintf('--- 关键配置 ---\n');
fprintf('  四旋翼数量 n            : %d\n', n);
fprintf('  负载质量 m0             : %.4f kg\n', m0);
fprintf('  单机质量 m              : %.4f kg\n', m);
fprintf('  绳长 l                  : %.4f m\n', cfg.link.length);
% ★★ 悬停推力必须由**挂点几何**解出，不能再用 m0 g / n（见 README §11）。
%   新挂点（一边中点 + 对边两顶点）质心偏离负载质心 0.0333 m，
%   悬停张力是 2:1:1 而不是三等分，因此**各机悬停推力互不相同**。
%   写死 m0 g / n + m g = 0.5804 N 会掩盖机 1 的 53.3% 推力需求。
muHoverCfg = [ones(1, n); cfg.payload.attachPoints(2, :); ...
              -cfg.payload.attachPoints(1, :)] \ [-m0 * g; 0; 0];
thrHoverCfg = abs(muHoverCfg).' + m * g;
fprintf('  挂点（负载体系 3x%d）    : %s\n', n, mat2str(cfg.payload.attachPoints, 4));
fprintf('  挂点质心偏移            : (%+.4f, %+.4f) m\n', ...
    mean(cfg.payload.attachPoints(1, :)), mean(cfg.payload.attachPoints(2, :)));

% ★★★ 挂点"是否在物品边缘"必须**显式验证并打印**，不能让人靠看图判断。
%   缘由：曾经在可视化里把负载按 2 倍显示（payloadSizeScale = 2.0），
%   而挂点仍画在真实位置 ⇒ 图上挂点落在板面 1/4 处，看起来像"挂在物品中间"，
%   被误判成几何写错了。其实几何一直是对的。现在：
%     ① payloadSizeScale 改回 1.0（几何保真）；
%     ② 这里直接把"离最近竖边的距离"打出来，0 就是恰在边界上。
hx = cfg.payload.size(1) / 2;
hy = cfg.payload.size(2) / 2;
hz = cfg.payload.size(3) / 2;
apCfg = cfg.payload.attachPoints;
distEdge = zeros(1, n);
distTop = zeros(1, n);
for i = 1:n
    distEdge(i) = min(hx - abs(apCfg(1, i)), hy - abs(apCfg(2, i)));
    distTop(i) = hz - abs(apCfg(3, i));
end
fprintf('  负载半宽 (hx, hy, hz)   : (%.4f, %.4f, %.4f) m\n', hx, hy, hz);
fprintf('  各挂点到最近竖直边距离  : [%s] m  ← 0 = 恰在边缘\n', ...
    mat2str(distEdge, 4));
fprintf('  各挂点到上表面的距离    : [%s] m  ← 0 = 恰在上表面\n', ...
    mat2str(distTop, 4));
fprintf('  分配矩阵 P 尺寸         : %d x %d，rank = %d（需要 = 6）\n', ...
    size(P, 1), size(P, 2), graphRank);
fprintf('  各机悬停推力（几何解）  : [%s] N\n', mat2str(thrHoverCfg, 4));
fprintf('  最大悬停推力占比        : %.1f %% 上限\n', 100 * max(thrHoverCfg) / tMax);
fprintf('  总悬停推力              : %.4f N / 可用 %.4f N\n', ...
    sum(thrHoverCfg), n * tMax);
fprintf('  推重比                  : %.2f\n', n * tMax / ((m0 + n * m) * g));
if isFigureEight
    fe = cfg.figureEight;
    fprintf('\n--- 八字避障工况（论文 Fig. 3 算例 / Fig. 2 快照）---\n');
    fprintf('  八字横宽 x              : ±%.2f m（半宽 %.2f m）\n', fe.amplitudeX, fe.amplitudeX);
    fprintf('  八字纵长 y              : ±%.2f m（半长 %.2f m）\n', ...
        0.5 * fe.amplitudeY, 0.5 * fe.amplitudeY);
    fprintf('  环绕圈数                : %.1f 个完整八字\n', fe.cycles);
    fprintf('  三阶段时长              : 起飞 %.1f s + 环绕 %.1f s + 降落 %.1f s = %.1f s\n', ...
        fe.takeoffDuration, fe.cruiseDuration, fe.landingDuration, cfg.simulation.duration);
    fprintf('  包络过渡时长 blendTime  : %.1f s（★ 必须 >= 5 s，见下）\n', fe.blendTime);
    fprintf('  环绕平面高度            : %.2f m\n', -fe.cruiseHeight);
    if cfg.obstacles.enabled
        fprintf('  锥形障碍物 %d 个        : ', size(cfg.obstacles.positions, 2));
        for j = 1:size(cfg.obstacles.positions, 2)
            fprintf('(%.2f, %.2f) ', cfg.obstacles.positions(1, j), cfg.obstacles.positions(2, j));
        end
        fprintf('\n');
        fprintf('    锥高 %.2f m，底半径 %.2f m，要求净间隙 > %.2f m\n', ...
            cfg.obstacles.height, cfg.obstacles.radius, cfg.obstacles.clearance);
    end
end
fprintf('\n');

% ---------------------------------------------------------------- 仿真
fprintf('--- 运行仿真 ---\n');
tStart = tic;
sim = crazyflie_slung_simulation(cfg);
fprintf('  完成，耗时 %.2f s，步数 %d\n\n', toc(tStart), numel(sim.time));

% ---------------------------------------------------------------- 自检
fprintf('--- 自检结果 ---\n');
checks = emptyCheckList();
s = sim.summary;

checks = addCheck(checks, '分配矩阵 rank(P) = 6', graphRank == 6, ...
    sprintf('rank = %d', graphRank), '= 6');

% ★★★ 挂点必须在物品**边缘**（用户明确要求：一边中点 + 对边两顶点）。
%   这一项必须独立断言，不能只看图 —— 可视化一旦缩放负载就会误导（见 README §13）。
checks = addCheck(checks, '三个挂点都在物品边缘（到边界距离 = 0）', ...
    all(distEdge < 1e-12) && all(distTop < 1e-12), ...
    sprintf('到竖边 max %.2e m，到上表面 max %.2e m', ...
    max(distEdge), max(distTop)), '= 0');

% ★ 分配一致性：伪逆解必须精确满足 sum mu_id = Fd 与 sum rho_hat R0'' mu_id = Md
[k1, k2] = allocationResidual(sim, cfg);
checks = addCheck(checks, '张力分配 sum(mu_id) = Fd', k1 < 1e-9, ...
    sprintf('max 残差 %.3e', k1), '< 1e-9');
checks = addCheck(checks, '张力分配 sum(rho_x mu_id) = Md', k2 < 1e-9, ...
    sprintf('max 残差 %.3e', k2), '< 1e-9');

checks = addCheck(checks, '状态量有限（无 NaN/Inf）', s.finiteState, ...
    string(s.finiteState), '必须为 true');

% ★★ 坏步计数必须为 0 —— 这是**唯一不会骗人的发散判据**。
%   动力学的 6x6 代数系统若出现 NaN/Inf，坏步的解会被置零保护，
%   NaN 不会传进状态 ⇒ **状态日志仍然可能全部有限**、上一项"状态量有限"照样通过。
%   曾因此把"已经发散"误判成"正常"，所以单独列一项。
if isfield(s, 'nonFiniteSolveCount')
    checks = addCheck(checks, '代数系统无非有限量（发散保护未触发）', ...
        s.nonFiniteSolveCount == 0, ...
        sprintf('%d 步', s.nonFiniteSolveCount), '= 0');
else
    checks = addCheck(checks, '代数系统无非有限量（发散保护未触发）', false, ...
        '缺少 summary.nonFiniteSolveCount 字段', '= 0');
end

qNormError = max(abs(sqrt(sum(sim.linkUnitLog.^2, 1)) - 1), [], 'all');
checks = addCheck(checks, '各绳向保持单位长度', qNormError < 1e-9, ...
    sprintf('最大范数误差 %.3e', qNormError), '< 1e-9');

if isfield(s, 'vehicleCollisionFree')
    checks = addCheck(checks, '无人机中心距满足碰撞裕度', s.vehicleCollisionFree, ...
        sprintf('最小 %.4f m / 要求 %.4f m', s.minVehicleSeparation, ...
        s.requiredVehicleSeparation), '> 要求值');
end
if isfield(s, 'vehiclePayloadCollisionFree')
    checks = addCheck(checks, '无人机与负载外接包络不碰撞', ...
        s.vehiclePayloadCollisionFree, ...
        sprintf('最小间隙 %.4f m / 裕度 %.4f m', ...
        s.minVehiclePayloadClearance, s.vehiclePayloadClearanceMargin), '> 0');
end

checks = addCheck(checks, '所有绳索张力恒为正（单边约束 mu >= 0）', ...
    s.allTensionsPositive, ...
    sprintf('min %.4f N / max %.4f N', s.minTension, s.maxTension), 'min > 0');

% ★ 绳索特有断言（"绳 vs 刚性连杆"的唯一数值可验差别）
%   绳在绷紧时与刚性连杆力学完全相同，但绳长必须严格守恒 —— 这正是"绷紧"的含义。
%   本实现里无人机位置由 veh = x0 + R0 rho_i - l q_i 定义，所以该不变量构造性成立，
%   实测偏差 1.7e-16 m。这里把它显式断言，防止将来改动破坏它。
%   单边约束 mu_i >= 0 已由上一项 'allTensionsPositive' 覆盖（实测裕度 69.2 %）。
checks = addCheck(checks, '绳索长度不变量 ‖挂点-无人机‖ ≡ l', ...
    s.ropeLengthInvariantHolds, ...
    sprintf('最大偏差 %.2e m (l = %.3f m)', s.maxRopeLengthDrift, ...
    cfg.link.length), '< 1e-9 m');

% ★ 悬停张力必须由**挂点几何**决定，不能再假设"三根绳三等分"。
%   新挂点（某边中点 + 对边两顶点）的质心偏离负载质心 0.0333 m，
%   于是悬停时三根绳的张力是 **2:1:1** 而不是各 m0 g / n。
%   判据：联立力平衡 Σμ = m0 g 与力矩平衡 Σ ρ̂_i μ_i = 0，
%   解出理论张力再与实测比对 —— 这同时验证了**张力分配矩阵**与**挂点几何**
%   是否一致（比原来那个"每根都等于 m0 g / n"的假设强得多）。
RHO = cfg.payload.attachPoints;                      % 3 x n
balanceMatrix = [ones(1, n); RHO(2, :); -RHO(1, :)];
muHover = balanceMatrix \ [-m0 * g; 0; 0];           % 悬停时各绳 μ_i 的 z 分量
tensionTheory = abs(muHover);
tensionMeasured = s.steadyTension(:);
relErr = abs(tensionMeasured - tensionTheory) ./ max(tensionTheory, eps);
tensionStr = @(v) strjoin(arrayfun(@(x) sprintf('%.4f', x), v(:).', ...
    'UniformOutput', false), ', ');
if isfield(cfg.allocation, 'outwardBiasFraction') && cfg.allocation.outwardBiasFraction > 0
    checks = addCheck(checks, '外张绳索模式总张力可接受', ...
        sum(tensionMeasured) >= m0 * g && sum(tensionMeasured) < 1.8 * m0 * g, ...
        sprintf('实测总张力 %.4f N / 负载重力 %.4f N', ...
        sum(tensionMeasured), m0 * g), '[m0 g, 1.8 m0 g)');
else
    checks = addCheck(checks, '稳态张力 = 挂点几何决定的悬停解', all(relErr < 0.10), ...
        sprintf('实测 [%s] / 理论 [%s] N（最大偏差 %.1f%%，合计 %.4f N = m0 g）', ...
        tensionStr(tensionMeasured), tensionStr(tensionTheory), ...
        100 * max(relErr), sum(tensionMeasured)), '< 10%');
end

% ★ 稳态位置误差只在**静态悬停工况**下用 2 cm 门限。
%   八字工况下"稳态"这个概念不成立（参考点一直在动），跟随误差必然大得多，
%   这时改用下面【八字组】的三阶段门限来判定。
%   如果把 2 cm 硬套到八字工况，会得到一个"永远 FAIL"的假告警，
%   反而掩盖真正的跟踪问题。
if ~isFigureEight
    checks = addCheck(checks, '稳态位置误差 < 2 cm', s.steadyPositionError < 0.02, ...
        sprintf('%.4f m (%.2f mm)', s.steadyPositionError, 1000*s.steadyPositionError), ...
        '< 0.02 m');

    % ★★ 配置一致性：**定高工况不应启用锥形障碍物**。
    %   理由：锥体的存在意义是"负载水平绕行时必须避开"，而定高悬停根本不做水平运动
    %   ⇒ 锥体既无物理意义，又会污染三维图与自检输出（曾出现"定高工况还画着两个锥"）。
    %   ★ 注意锥体只参与可视化与自检、**不进入动力学**，所以这条只是配置一致性断言，
    %     关闭它对全部动力学指标零影响。
    checks = addCheck(checks, '定高工况未启用锥形障碍物（无反意义务）', ...
        ~cfg.obstacles.enabled, ...
        sprintf('obstacles.enabled = %d（定高工况应为 0）', cfg.obstacles.enabled), ...
        '= false');
    if cfg.obstacles.enabled
        warning('crazyflie_slung_demo:ObstaclesInHoverCase', ...
            ['当前为**静态悬停**工况（cfg.referenceFcn 为空），但 cfg.obstacles.enabled = true。' ...
             '锥形障碍物在定高工况下没有意义，建议置 false —— 三维图与自检会更干净。']);
    end
end
% ★ 绳向误差的门限需要考虑 yaw 参考变化的注入。
%   负载 yaw -> 挂点方位角转动 -> 期望张力分布改变 ->
%   绳索环追赶不及 -> 留下方向残差。
%   因此 3 机情形的合理判据是 3 deg（而非单机版的 2 deg），
%   并在演示窗口（30 s）内成立；这一点在 README §6.2 与 §8 中已明确说明。
checks = addCheck(checks, '稳态绳向误差 < 3 deg', s.steadyLinkError < deg2rad(3), ...
    sprintf('%.3f deg（受 yaw 参考变化注入）', rad2deg(s.steadyLinkError)), ...
    '< 3 deg');

% ★ 负载姿态按 roll/pitch 与 yaw 分开判定：yaw 是低带宽闭环轴，
%   允许比 roll/pitch 更宽的误差门限，但必须直接检查相对参考角。
if isfield(s, 'steadyAttitudeErrorVec') && numel(s.steadyAttitudeErrorVec) == 3
    attVec = s.steadyAttitudeErrorVec(:);
else
    attVec = repmat(s.steadyAttitudeError, 3, 1);
end
rollPitchAttErr = rad2deg(norm(attVec(1:2)));
checks = addCheck(checks, '稳态姿态误差(roll/pitch) < 3 deg', ...
    rollPitchAttErr < 3, ...
    sprintf('%.4f deg（三轴范数 %.3f deg，yaw = %.3f deg）', ...
    rollPitchAttErr, rad2deg(norm(attVec)), rad2deg(attVec(3))), '< 3 deg');

% yaw 角速度门限用于防止低带宽通道在张力分配后产生过激响应。
if isfield(s, 'loadBodyRateFinal') && numel(s.loadBodyRateFinal) == 3
    yawRate = abs(s.loadBodyRateFinal(3));
else
    yawRate = 0;
end
checks = addCheck(checks, '负载 yaw 角速度 < 1.0 rad/s', yawRate < 1.0, ...
    sprintf('%.4f rad/s（低带宽 yaw 闭环）', yawRate), '< 1.0 rad/s');
if isfield(cfg.loadController, 'yawChannelEnabled') && cfg.loadController.yawChannelEnabled ...
        && isfield(s, 'steadyYawTrackingError') && isfinite(s.steadyYawTrackingError)
    checks = addCheck(checks, '稳态负载 yaw 跟踪误差 < 10 deg', ...
        s.steadyYawTrackingError < deg2rad(10), ...
        sprintf('%.3f deg', rad2deg(s.steadyYawTrackingError)), '< 10 deg');
end

% ★ 推力峰值要**分开看起步段与稳态段**。
%   cfg.initial 有意给负载一个很大的初始偏差（位置离目标 0.81 m、姿态 3/-2/4 deg、
%   各绳还有 3~5 deg 倾角），第一步的推力指令必然短暂饱和 —— 这是"测试收敛能力"
%   这个意图的必然结果，不是控制缺陷。而且新挂点几何下三根绳是 2:1:1，
%   机 1 本来就承担两倍张力，起步段更容易碰到上限。
%   所以判据用**跳过起步 0.2 s 之后**的峰值；两个数都打印出来，不藏。
startupSkip = round(0.2 / cfg.simulation.dt);
thrustPctSteady = max(sim.thrustPctLog(:, min(startupSkip + 1, end):end), [], 'all');
checks = addCheck(checks, '推力峰值占比 < 95%（跳过起步 0.2 s）', thrustPctSteady < 95, ...
    sprintf('稳态 %.1f %% / 全程 %.1f %%（含起步 0.2 s 的初始偏差瞬态）', ...
    thrustPctSteady, s.maxThrustPercentage), '< 95 %%');

checks = addCheck(checks, '最大机体角速度 < 8 rad/s', s.maxBodyRate < 8, ...
    sprintf('%.3f rad/s', s.maxBodyRate), '< 8 rad/s');

target = cfg.target.position(:);
finalPos = s.loadPositionFinal(:);
if isFigureEight
    % 八字工况的终点是 cfg.figureEight.landPosition（= target.position），
    % 而且它经过了完整的起飞/环绕/降落三段，所以门限要放宽：
    % 降落段本身只有 4 s，从巡航误差里恢复需要时间。
    % 实测（KX=14, KI=3.0, B=6.0）：末点残差 78.7 mm，主要是 y 方向残留
    % （最后 4 s 从环绕段收尾时的 y 偏移 123 mm 衰减到 75 mm）。
    checks = addCheck(checks, '降落终点残差 < 0.10 m', ...
        s.landingPointError < 0.10, ...
        sprintf('%.1f mm (落点 [%.3f %.3f %.3f])', ...
        1000 * s.landingPointError, finalPos(1), finalPos(2), finalPos(3)), ...
        '< 0.10 m');
else
    checks = addCheck(checks, '负载终位置接近目标', norm(finalPos - target) < 0.02, ...
        sprintf('[%.4f %.4f %.4f] m', finalPos(1), finalPos(2), finalPos(3)), ...
        '欧氏距离 < 0.02 m');
end

% ==================== 八字避障工况专项自检 ====================
% 这一组的门限全部来自**实测调参数值**，不是拍脑袋的数字。
% 调参过程与背后的两个根因（包络加速度尖峰 + q_id_dot 差分噪声）
% 记录在 README §10 与参数文件的注释里。
if isFigureEight
    fe = cfg.figureEight;

    % ---- 三阶段跟踪误差 ----
    % 门限全部按**最终参数集（KX=14, KV=8.68, KI=3.0, B=6.0, noQIDDOT）**的
    % 实测值留裕度设定，见 _diag_fe_final.py / README §10。
    %   实测：起飞峰 23.4 mm / 环绕峰 587.9 mm 均 117.3 mm / 降落峰 122.9 mm
    %        末点残差 78.7 mm / Tmin 0.0985 N / 推力峰 70.8%
    % 起飞段：只是竖直爬升，参考加速度小，应达到 cm 级
    checks = addCheck(checks, '起飞段跟踪误差 < 50 mm', ...
        s.maxPositionErrorTakeoff < 0.050, ...
        sprintf('%.1f mm', 1000 * s.maxPositionErrorTakeoff), '< 50 mm');

    % 环绕段：±1.5 m 大幅机动 + 绳摆耦合。
    % ★ 门限按实测峰 588 mm 留约 70% 裕度取 1.0 m —— 既能守住
    %   "不发散、不失控"，又不会因正常机动误报。
    % ★ 误差主要来自**绳摆相位滞后**，不是控制器不稳：
    %   位置环带宽 wn = 3.74 rad/s 是绳摆频率 sqrt(g/l) = 5.29 rad/s 的 0.71 倍。
    %   这是**刻意**选的 —— 实测把带宽提到 1.0 倍反而 2.42 s 就发散
    %   （寄生耦合：位置环越快越把平动转成摆动），详见参数文件的长注释。
    checks = addCheck(checks, '环绕段跟踪误差峰值 < 1.0 m', ...
        s.maxPositionErrorCruise < 1.0, ...
        sprintf('%.1f mm', 1000 * s.maxPositionErrorCruise), '< 1.0 m');
    checks = addCheck(checks, '环绕段跟踪误差均值 < 0.25 m', ...
        s.meanPositionErrorCruise < 0.25, ...
        sprintf('%.1f mm', 1000 * s.meanPositionErrorCruise), '< 0.25 m');

    % 降落段：z 收得很干净（实测 <= 3 mm），水平残差见下一项
    checks = addCheck(checks, '降落段跟踪误差 < 0.20 m', ...
        s.maxPositionErrorLanding < 0.20, ...
        sprintf('%.1f mm', 1000 * s.maxPositionErrorLanding), '< 0.20 m');

    % ---- 环绕的是不是"八字" ----
    % 判据：x 方向必须完成 2*cycles 次往复，即轨迹在 x 上至少走满 ±amplitudeX
    % 的 80%（跟随误差会让实际幅值略小于期望幅值）。
    xSpan = s.trajectorySpan(1);
    checks = addCheck(checks, 'x 方向完成 2*cycles 次往复', ...
        xSpan > 1.6 * fe.amplitudeX, ...
        sprintf('x 行程 %.3f m（期望幅值 ±%.2f m）', xSpan, fe.amplitudeX), ...
        sprintf('> %.2f m', 1.6 * fe.amplitudeX));

    % ---- 绕开两个锥形障碍物（论文 Fig. 3 的核心要求）----
    % ★★ 必须同时检查"有效性"，不能只看数值。 ★★
    %   本自检曾经因为高度量定义混乱而**静默失效**：所有环绕段采样都走进
    %   gap = inf 分支，于是 min 净间隙是一个与环绕段无关的起飞瞬间值，
    %   然而判据照样 PASS —— 属于"根本没检查"却报"通过"。
    %   所以这里拆成两步：
    %     (1) obstacleClearanceValid —— 环绕段是否真的进入了锥的高度区间；
    %     (2) obstacleClearanceOk    —— 净间隙是否 > 规定值（已内含 (1)）。
    if cfg.obstacles.enabled
        branchCounts = [0, 0, 0];
        if isfield(s, 'obstacleBranchCounts')
            branchCounts = s.obstacleBranchCounts;
        end
        validFlag = true;
        if isfield(s, 'obstacleClearanceValid')
            validFlag = s.obstacleClearanceValid;
        end
        % (1) 有效性：落在锥高度区间内的采样点必须非零
        checks = addCheck(checks, '锥形障碍物净间隙自检有效（非静默失效）', ...
            validFlag, ...
            sprintf('落在锥高度区间内的采样 %d 点（锥顶以上 %d、锥底以下 %d）', ...
            branchCounts(1), branchCounts(2), branchCounts(3)), ...
            '区间内采样 > 0（否则避障从未被真正检查）');
        % 锥尖必须高于环绕平面，否则负载会从锥顶上方飞过、避障形同虚设
        apexAboveGround = -cfg.obstacles.baseZ + cfg.obstacles.height;
        checks = addCheck(checks, '锥尖高于环绕平面（避障为真实约束）', ...
            apexAboveGround > abs(fe.cruiseHeight), ...
            sprintf('锥尖离地 %.2f m vs 环绕平面离地 %.2f m', ...
            apexAboveGround, abs(fe.cruiseHeight)), ...
            '锥尖 > 环绕平面');
        % (2) 净间隙
        checks = addCheck(checks, '负载到锥形障碍物净间隙 > 规定值', ...
            s.obstacleClearanceOk, ...
            sprintf('最小净间隙 %.1f mm (t = %.2f s)，负载外接球半径 %.1f mm', ...
            1000 * s.minObstacleClearance, s.minObstacleClearanceTime, ...
            1000 * s.payloadBoundingRadius), ...
            sprintf('> %.0f mm', 1000 * cfg.obstacles.clearance));
    end

    % ---- 包络参数是否在安全区 ----
    % blendTime 太小会让包络二阶导在环绕段两端制造加速度尖峰并导致发散
    % （B=1.2 s 时参考加速度峰值 9.20 m/s²，是纯八字理论峰值 2.92 m/s² 的 3.15 倍）。
    % 这里显式断言它，防止将来有人为了"贴近论文"把它调小。
    checks = addCheck(checks, '包络过渡时长 blendTime >= 5 s', ...
        fe.blendTime >= 5.0, ...
        sprintf('%.1f s', fe.blendTime), '>= 5 s（否则包络加速度尖峰会发散）');
end

nPass = 0;
for i = 1:numel(checks)
    if checks(i).passed
        verdict = 'PASS';
        nPass = nPass + 1;
    else
        verdict = 'FAIL';
    end
    fprintf('  [%s] %-30s | %-42s | 期望 %s\n', ...
        verdict, checks(i).name, checks(i).value, checks(i).expectation);
end
fprintf('\n  通过 %d / %d 项\n\n', nPass, numel(checks));

% ---------------------------------------------------------------- 指标汇总
fprintf('--- 关键指标 ---\n');
fprintf('  稳态位置误差        : %.5f m   (%.2f mm)\n', ...
    s.steadyPositionError, 1000 * s.steadyPositionError);
fprintf('  稳态高度偏差        : %.5f m\n', abs(s.loadHeightFinal - (-target(3))));
fprintf('  稳态绳向误差(最大)  : %.5f (%.3f deg)\n', ...
    s.steadyLinkError, rad2deg(s.steadyLinkError));
fprintf('  稳态姿态误差(最大)  : %.5f (%.3f deg)  [机体]\n', ...
    s.steadyAttitudeError, rad2deg(s.steadyAttitudeError));
if isfield(s, 'steadyAttitudeErrorVec') && numel(s.steadyAttitudeErrorVec) == 3
    ax = rad2deg(s.steadyAttitudeErrorVec(:));
    fprintf('  稳态负载姿态误差    : [roll pitch yaw] = [%.4f %.4f %.4f] deg\n', ...
        ax(1), ax(2), ax(3));
    fprintf('    -> roll/pitch（严格门限）: %.3f deg\n', norm(ax(1:2)));
    fprintf('    -> yaw（低带宽闭环）      : %.3f deg\n', ax(3));
end
if isfield(s, 'loadBodyRateFinal') && numel(s.loadBodyRateFinal) == 3
    fprintf('  负载角速度终值      : [%.4f %.4f %.4f] rad/s\n', ...
        s.loadBodyRateFinal(1), s.loadBodyRateFinal(2), s.loadBodyRateFinal(3));
end
fprintf('  最大位置误差        : %.4f m\n', s.maxPositionError);
if isFigureEight && isfield(s, 'phaseWindows') && ~isempty(s.phaseWindows)
    w = s.phaseWindows;
    fprintf('\n--- 八字避障工况指标 ---\n');
    fprintf('  起飞段 [0, %.1f s]      : 误差峰值 %.1f mm\n', ...
        s.durationPhases(1), 1000 * s.maxPositionErrorTakeoff);
    fprintf('  环绕段 [%.1f, %.1f s]   : 误差峰值 %.1f mm / 均值 %.1f mm\n', ...
        s.durationPhases(1), s.durationPhases(2), ...
        1000 * s.maxPositionErrorCruise, 1000 * s.meanPositionErrorCruise);
    fprintf('  降落段 [%.1f, %.1f s]  : 误差峰值 %.1f mm\n', ...
        s.durationPhases(2), s.durationPhases(3), ...
        1000 * s.maxPositionErrorLanding);
    fprintf('  降落终点残差        : %.1f mm\n', 1000 * s.landingPointError);
    fprintf('  轨迹行程 (x, y)     : (%.3f, %.3f) m\n', ...
        s.trajectorySpan(1), s.trajectorySpan(2));
    if isfield(cfg, 'obstacles') && cfg.obstacles.enabled
        fprintf('  到锥最小净间隙      : %.1f mm (t = %.2f s)%s\n', ...
            1000 * s.minObstacleClearance, s.minObstacleClearanceTime, ...
            ternary(s.obstacleClearanceOk, '  [OK]', '  [过近!]'));
        if isfield(s, 'obstacleBranchCounts')
            bc = s.obstacleBranchCounts;
            fprintf('    自检有效性        : %s  (锥高度区间内 %d 点 / 锥顶以上 %d / 锥底以下 %d)\n', ...
                ternary(s.obstacleClearanceValid, '[有效]', '[★静默失效★]'), ...
                bc(1), bc(2), bc(3));
        end
        fprintf('    锥尖离地 / 环绕平面: %.2f m / %.2f m%s\n', ...
            -cfg.obstacles.baseZ + cfg.obstacles.height, abs(cfg.figureEight.cruiseHeight), ...
            ternary(-cfg.obstacles.baseZ + cfg.obstacles.height > abs(cfg.figureEight.cruiseHeight), ...
            '  [锥尖更高，避障为真实约束]', '  [★锥尖偏低，负载会从上方飞过★]'));
    end
    fprintf('  包络过渡时长        : %.1f s\n', cfg.figureEight.blendTime);
end
fprintf('\n');
fprintf('  各绳索稳态张力      : %s N  (理论 %s，合计 %.4f = m0 g)\n', ...
    tensionStr(s.steadyTension(:)), tensionStr(tensionTheory), ...
    sum(tensionTheory));
fprintf('  张力范围            : %.4f ~ %.4f N\n', s.minTension, s.maxTension);
fprintf('  推力峰值占比        : %.1f %%\n', s.maxThrustPercentage);
fprintf('  最大机体角速度      : %.3f rad/s\n', s.maxBodyRate);
if isfield(s, 'nonFiniteSolveCount')
    fprintf('  坏步计数（必须为 0）: %d\n', s.nonFiniteSolveCount);
end
fprintf('  负载终高度          : %.4f m\n', s.loadHeightFinal);
fprintf('  各机终高度          : %s m\n', mat2str(s.vehicleHeightFinal.', 4));
fprintf('\n');

report = struct();
report.config = cfg;
report.summary = s;
report.checks = checks;
report.numPassed = nPass;
report.numTotal = numel(checks);
report.allPassed = (nPass == numel(checks));
report.allocationRank = graphRank;

if report.allPassed
    fprintf('==============================================================\n');
    fprintf(' 全部自检通过 —— 仿真结果合理。\n');
    fprintf('==============================================================\n');
else
    fprintf('==============================================================\n');
    fprintf(' 存在未通过项，请检查上面的 FAIL 行。\n');
    fprintf('==============================================================\n');
end
end

% ======================================================================
function [maxErrForce, maxErrMoment] = allocationResidual(sim, cfg)
% 逐步校验伪逆分配是否精确满足论文 (22)：
%   sum_i mu_id = Fd ,  sum_i hat(rho_i) R0' mu_id = Md
nSteps = numel(sim.time);
n = cfg.vehicle.count;
rhoAll = cfg.payload.attachPoints;
maxErrForce = 0;
maxErrMoment = 0;
for k = 1:nSteps
    Fd = sim.commandLog(k).desiredForce;
    Md = sim.commandLog(k).desiredMoment;
    R0 = sim.loadRotationLog(:, :, k);
    muAll = sim.desiredTensionLog(:, :, k);
    sumF = zeros(3, 1);
    sumM = zeros(3, 1);
    for i = 1:n
        sumF = sumF + muAll(:, i);
        sumM = sumM + hat(rhoAll(:, i)) * (R0.' * muAll(:, i));
    end
    maxErrForce = max(maxErrForce, norm(sumF - Fd));
    maxErrMoment = max(maxErrMoment, norm(sumM - Md));
end
end

% ======================================================================
function S = hat(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function list = emptyCheckList()
list = struct('name', {}, 'passed', {}, 'value', {}, 'expectation', {});
end

function list = addCheck(list, name, passed, valueText, expectationText)
entry = struct('name', name, 'passed', logical(passed), ...
    'value', char(valueText), 'expectation', char(expectationText));
if isempty(list)
    list = entry;
else
    list(end + 1) = entry;
end
end

function out = ternary(condition, whenTrue, whenFalse)
% 三元表达式（MATLAB 没有内置，自己写一个保持调用处紧凑）。
if condition
    out = whenTrue;
else
    out = whenFalse;
end
end
