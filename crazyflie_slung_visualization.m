function crazyflie_slung_visualization(sim, cfg)
%CRAZYFLIE_SLUNG_VISUALIZATION 多机吊运仿真的结果绘图与三维动画。
%
% 本文件不参与任何控制或动力学计算，只消费 sim 中的日志数据。

if cfg.visualization.plot
    plotSummary(sim, cfg);
end

if cfg.visualization.animate
    animateScene(sim, cfg);
end
end

% ======================================================================
function plotSummary(sim, cfg)
% 汇总曲线图（6 个子图）。
time = sim.time;
nSteps = numel(time);
n = cfg.vehicle.count;
steady = max(1, round(0.8 * nSteps)):nSteps;

% 显示时把 z 取反（绘图纵轴向上为正）
loadDisplay = sim.loadPositionLog;
loadDisplay(3, :) = -loadDisplay(3, :);
vehDisplay = squeeze(sim.vehiclePositionLog);      % 3 x n x nSteps
if n == 1
    vehDisplay = reshape(vehDisplay, 3, 1, nSteps);
end
vehDisplay(3, :, :) = -vehDisplay(3, :, :);

figure('Name', '多机协同吊运仿真结果', 'Color', 'w', ...
    'Position', [70, 50, 1250, 740]);

% ---- (1) 三维轨迹 ----
% ★ 轴框**不能**按负载轨迹自动适应。原因：起飞时负载从地面接触高度爬到
%   运输高度；如果把初始起飞段和异常过冲一起用于定框，负载边长和绳索都会
%   被压缩成很小的视觉对象，无法判断是否穿地或是否进入绷紧阶段。
%   改为以**四旋翼集群在稳态窗口的位置**为中心定框，负载爬升的那一段
%   会自然伸出框外（视觉上仍然完整可读，因为它是单调上升的一条线）。
%   子图只画前 trajectoryWindow 秒，避免起飞/收紧瞬态占满画面。
subplot(2, 3, 1);
if isfield(cfg, 'takeoff') && cfg.takeoff.enabled ...
        && cfg.takeoff.landingEnabled
    % 起降混合模式下必须把末段 LANDING_RELEASE 也纳入三维轨迹，
    % 否则用户只能看到起飞和运输，看不到无人机解除绳索后的降落。
    winEnd = nSteps;
else
    winEnd = min(nSteps, max(2, round(cfg.visualization.trajectoryWindow / ...
        cfg.simulation.dt) + 1));
end
winSel = 1:winEnd;
plot3(loadDisplay(1, winSel), loadDisplay(2, winSel), loadDisplay(3, winSel), ...
    'LineWidth', 1.8, 'Color', [0.00, 0.35, 0.75]);
hold on; grid on; axis equal;
vehColors = lineColors(n);
for i = 1:n
    plot3(squeeze(vehDisplay(1, i, winSel)), squeeze(vehDisplay(2, i, winSel)), ...
        squeeze(vehDisplay(3, i, winSel)), 'LineWidth', 1.0, 'Color', vehColors(i, :));
end
plot3(loadDisplay(1, winSel(end)), loadDisplay(2, winSel(end)), ...
    loadDisplay(3, winSel(end)), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.00, 0.35, 0.75], 'Color', 'none');

% 轴框：以**稳态窗口内四旋翼集群的质心**为中心（用截面点，不用整段轨迹）
%
% ★★ 数值安全（必须保留）：若仿真中途发散，状态里会出现 NaN/Inf。
%    此时下面的"稳态窗口质心 + 跨度"会算成 NaN，接着
%    xlim/ylim/zlim 会直接抛
%        "范围必须为包含递增的数值的 2 元素向量"
%    把整个结果展示打断 —— 用户就看不到任何指标，只看到一个报错。
%    这里改成：只统计**有限样本**；若稳态窗口全非有限，退回整段的有限样本；
%    若整段都没有有限样本，退回默认轴框。曲线里的断点即发散时刻，仍然可见。
steadySel = max(1, round(0.8 * nSteps)):nSteps;
vehCluster = squeeze(mean(vehDisplay(:, :, steadySel), 2));    % 3 x nSteady
okCols = all(isfinite(vehCluster), 1);
if ~any(okCols)
    stepFinite = reshape(squeeze(all(all(isfinite(vehDisplay), 1), 2)), 1, []);
    if numel(stepFinite) == nSteps && any(stepFinite)
        vehCluster = squeeze(mean(vehDisplay(:, :, stepFinite), 2));
        okCols = all(isfinite(vehCluster), 1);
    end
end
if any(okCols)
    vehCluster = vehCluster(:, okCols);
    clusterMid = mean(vehCluster, 2);
    % 跨度取"四旋翼集群自身尺度"与最小跨度中的较大者，保证机臂/绳索看得清
    clusterSpan = max([max(vehCluster(1, :)) - min(vehCluster(1, :)), ...
        max(vehCluster(2, :)) - min(vehCluster(2, :))]);
else
    % 全段都没有有限样本 ⇒ 退回默认框（画面上会是一条空白 + 断开曲线）
    warning('crazyflie_slung:NonFiniteState', ...
        ['轨迹中没有任何有限的状态样本（仿真从很早就不发散地失败了），' ...
         '三维轴框退回默认范围。']);
    clusterMid = [0; 0; 0];
    clusterSpan = cfg.visualization.trajectoryMinSpan;
end
trjHalf = max([cfg.visualization.trajectoryMinSpan, ...
    2.0 * clusterSpan, 2.2 * norm(cfg.payload.size)]) / 2;
if ~isfinite(trjHalf) || trjHalf <= 0
    trjHalf = cfg.visualization.trajectoryMinSpan / 2;
end
xlim([clusterMid(1) - trjHalf, clusterMid(1) + trjHalf]);
ylim([clusterMid(2) - trjHalf, clusterMid(2) + trjHalf]);
zlim([clusterMid(3) - trjHalf, clusterMid(3) + trjHalf]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('高度 (m)');
title(sprintf('负载与 %d 架四旋翼轨迹（%.0f s，稳态局部放大）', n, ...
    time(winEnd)));
legendEntries = [{'负载'}, arrayfun(@(i) sprintf('四旋翼 %d', i), 1:n, ...
    'UniformOutput', false), {sprintf('t=%.0fs 处', time(winEnd))}];
legend(legendEntries{:}, 'Location', 'best');

% ---- (2) 负载位置误差 ----
subplot(2, 3, 2);
plot(time, sim.positionErrorVectorLog.', 'LineWidth', 1.2);
grid on; xlabel('时间 (s)'); ylabel('误差 (m)');
title('负载位置误差');
legend('e_x', 'e_y', 'e_z', 'Location', 'best');

% ---- (3) 各绳向误差范数 ----
subplot(2, 3, 3);
plot(time, rad2deg(sim.linkErrorLog).', 'LineWidth', 1.2);
grid on; xlabel('时间 (s)'); ylabel('‖e_q‖ (deg)');
title('各绳向误差 ‖q_{id} × q_i‖');
legend(arrayfun(@(i) sprintf('绳索 %d', i), 1:n, 'UniformOutput', false), ...
    'Location', 'best');

% 注意：该误差**没有真正的稳态值**。它会随负载偏航漂移缓慢单调增长
% （实测斜率 +0.057 deg/s），自检里的 3 deg 门限对应 30 s 的仿真时长。
% 详见 README §6.2。

% ---- (4) 高度 ----
subplot(2, 3, 4);
plot(time, -loadDisplay(3, :), 'LineWidth', 1.6, 'Color', [0.00, 0.35, 0.75]);
hold on; grid on;
for i = 1:n
    plot(time, squeeze(vehDisplay(3, i, :)), 'LineWidth', 1.0, ...
        'Color', vehColors(i, :));
end
plot([time(1), time(end)], repmat(-cfg.target.position(3), 1, 2), '--', ...
    'Color', [0.85, 0.25, 0.10], 'LineWidth', 1.0);
plot([time(1), time(end)], [0, 0], ':', ...
    'Color', [0.15, 0.15, 0.15], 'LineWidth', 1.0);
xlabel('时间 (s)'); ylabel('高度 (m)');
title('高度（负载 vs 各四旋翼）');
legend([{'负载'}, arrayfun(@(i) sprintf('机 %d', i), 1:n, 'UniformOutput', false), ...
    {'目标高度', '地面'}], 'Location', 'best');

% ---- (5) 各绳索张力（左轴）与各机推力（右轴）----
subplot(2, 3, 5);
yyaxis left;
plot(time, sim.tensionLog.', 'LineWidth', 1.2);
hold on;
plot([time(1), time(end)], ...
    repmat(cfg.payload.mass * cfg.vehicle.gravity / n, 1, 2), ':', ...
    'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0);
ylabel('绳索张力 (N)');
yyaxis right;
hThrust = plot(time, sim.thrustPctLog.', '--', 'LineWidth', 1.0);
ylabel('推力占比 (%)');
hold on;
plot([time(1), time(end)], repmat(100, 1, 2), ':', 'Color', [0.6, 0.6, 0.6]);
grid on; xlabel('时间 (s)');
title('绳索张力（左）与推力占比（右）');
legend(hThrust, arrayfun(@(i) sprintf('机 %d 推力', i), 1:n, ...
    'UniformOutput', false), 'Location', 'best');

% ---- (6) 负载偏航角（yaw 跟踪通道）----
% ★ 这一格原来画"机 1 角速度指令 vs 实测"，但它与自检输出的
%   "最大机体角速度" 完全重复。负载偏航角才是本仿真最需要被看见的量：
%   它直接反映负载 yaw 参考与实际姿态的闭环跟踪结果。
subplot(2, 3, 6);
% ★ 同时画【参考 yaw】与【实际 yaw】。参考取仿真时记下的 sim.loadYawRefLog
%   （= 参考姿态 R0d 第一轴方位角：八字 = 运动方向；定高恒为 0）。
%   ★★ 两条都必须 unwrap：yaw 是角度，+179° 跳到 -179° 只是跨过 ±180 割线、
%      并非真的反向；不解卷绕就会被看成"正负乱跳"。
%   ★★ 有效性判据只看**字段是否存在、长度是否对得上**；
%      **绝不能用 any(loadYawRefLog ~= 0)** —— 定高工况参考 yaw 恒为 0，
%      用后者会把"恒为 0 的合法数据"误判成"无数据"，参考曲线就不画了。
yawDeg = rad2deg(atan2(squeeze(sim.loadRotationLog(2, 1, :)), ...
    squeeze(sim.loadRotationLog(1, 1, :))));
yawUnwrap = unwrapLocal(yawDeg);

yawDesUnwrap = [];
if isfield(sim, 'loadYawRefLog') && numel(sim.loadYawRefLog) == numel(time)
    yawDesUnwrap = unwrapLocal(rad2deg(sim.loadYawRefLog(:).'));
    % 两条曲线对齐到同一支：各自 unwrap 后可能相差整数个 360°，
    % 不对齐会出现"一条 +170、另一条 −190"的假象。
    yawUnwrap = yawUnwrap + 360 * round((yawDesUnwrap(1) - yawUnwrap(1)) / 360);
end

hold on; grid on;
if ~isempty(yawDesUnwrap)
    plot(time, yawDesUnwrap, '--', 'LineWidth', 1.8, 'Color', [0.15, 0.45, 0.75]);
end
plot(time, yawUnwrap, '-', 'LineWidth', 1.6, 'Color', [0.55, 0.25, 0.65]);
xlabel('时间 (s)'); ylabel('偏航角 (deg，已解卷绕)');
if ~isempty(yawDesUnwrap)
    title(sprintf('负载偏航角：参考(虚线) vs 实际(实线)  |  跟踪误差 %.1f°rms', ...
        sqrt(mean((yawUnwrap - yawDesUnwrap).^2))));
    legend({'参考 yaw', '实际 yaw'}, 'Location', 'best');
else
    title('负载偏航角（yaw 反馈关闭，仅作观测）');
    legend({'实际 yaw'}, 'Location', 'best');
end

% 起飞/收紧/降落诊断单独成图，避免把绳索松弛时的零张力误读为控制失败。
if isfield(sim, 'takeoffModeLog') && isfield(sim, 'ropeSlackLog')
    figure('Name', '绳索状态与起降阶段', 'Color', 'w', ...
        'Position', [160, 120, 980, 430]);
    subplot(2, 1, 1);
    plot(time, sim.ropeSlackLog.', 'LineWidth', 1.1);
    hold on; grid on;
    plot(time, zeros(size(time)), 'k:', 'LineWidth', 1.0);
    ylabel('绳长余量 l-d_i (m)');
    title('实际定位几何估计的绳索余量（正值=松弛，0=绷紧）');
    legend(arrayfun(@(i) sprintf('绳索 %d', i), 1:n, 'UniformOutput', false), ...
        'Location', 'best');
    subplot(2, 1, 2);
    stairs(time, sim.takeoffModeLog, 'LineWidth', 1.3, ...
        'Color', [0.20, 0.35, 0.70]);
    hold on; grid on;
    plot(time, sim.tensionScaleLog, '--', 'Color', [0.85, 0.25, 0.15], ...
        'LineWidth', 1.1);
    xlabel('时间 (s)'); ylabel('阶段编码 / 张力比例');
    yticks(0:5);
    yticklabels({'SLACK', 'TAKEUP', 'RAMP', 'ACTIVE', 'LAND-Taut', 'LAND-Release'});
    ylim([-0.3, 5.3]);
    legend({'阶段', '张力软启动比例'}, 'Location', 'best');
end

% 各机**机体**姿态角单独成图。
% ★ 与上面第 6 格的"负载偏航角"不是一回事：那个是**负载**的姿态，
%   这里是每架四旋翼**自己**的 roll/pitch/yaw。
plotVehicleAttitude(sim, cfg);
end

% ======================================================================
function plotVehicleAttitude(sim, cfg)
%PLOTVEHICLEATTITUDE 各架四旋翼的机体姿态角（roll / pitch / yaw）随时间变化。
%
% 三个子图各一条时间曲线×n 架（本工况 n=3 ⇒ 共 3x3 = 9 条）。
%
% 姿态角由 rotationLog 中的旋转矩阵按 **ZYX 欧拉角**反解，与
% crazyflie_slung_parameters.m 的 rpyToRotm 是同一个约定
% （R = Rz(yaw) * Ry(pitch) * Rx(roll)）：
%     pitch = atan2( -R(3,1), sqrt(R(1,1)^2 + R(2,1)^2) )
%     roll  = atan2(  R(3,2), R(3,3) )
%     yaw   = atan2(  R(2,1), R(1,1) )
% 这里用 atan2 而不是 asin(-R(3,1))，避免俯仰接近 ±90° 时的数值退化。
%
% yaw 在 ±180° 处会跳变，做一次解卷绕（复用本文件的 unwrapLocal）；
% roll / pitch 不跨割线，不需要解卷绕。

if ~isfield(sim, 'rotationLog')
    return      % 旧日志没有该字段时安静跳过，不影响其它图
end

n = cfg.vehicle.count;
time = sim.time;
nSteps = numel(time);
R = sim.rotationLog;                       % 3 x 3 x n x nSteps

% ★ 三轴各自存成 **n x nSteps 的二维矩阵**，不要写成 `rpy(1, :, :) = <n x nSteps>`：
%   那种"把二维矩阵赋给三维切片"的写法形状不匹配（1xnxnSteps vs nxnSteps），
%   MATLAB 的行为依赖版本/上下文，容易变成静默广播。全用二维就没有歧义。
%   （R(3,2,:,:) 是 [1 1 n nSteps]，按列主序 reshape 成 [n nSteps] 后
%     行 = 机体序号、列 = 时间步，已用 Python 逐元素验证过。）
rollSeries  = rad2deg(atan2(reshape(R(3, 2, :, :), n, nSteps), ...
                            reshape(R(3, 3, :, :), n, nSteps)));
pitchSeries = rad2deg(atan2(-reshape(R(3, 1, :, :), n, nSteps), ...
    sqrt(reshape(R(1, 1, :, :), n, nSteps).^2 ...
       + reshape(R(2, 1, :, :), n, nSteps).^2)));
yawSeries   = rad2deg(atan2(reshape(R(2, 1, :, :), n, nSteps), ...
                            reshape(R(1, 1, :, :), n, nSteps)));
for i = 1:n
    yawSeries(i, :) = unwrapLocal(yawSeries(i, :));
end
series = {rollSeries, pitchSeries, yawSeries};

titles = {'滚转角 roll', '俯仰角 pitch', '偏航角 yaw（已解卷绕，仅作观测）'};
colors = lineColors(n);
figure('Name', '四旋翼机体姿态角（roll / pitch / yaw）', 'Color', 'w', ...
    'Position', [200, 60, 900, 760]);
for a = 1:3
    subplot(3, 1, a);
    hold on; grid on;
    for i = 1:n
        plot(time, series{a}(i, :), 'LineWidth', 1.1, 'Color', colors(i, :));
    end
    ylabel(sprintf('%s (°)', titles{a}));
    title(sprintf('%s   |   %d 架中最大绝对值 %.2f°', ...
        titles{a}, n, max(abs(series{a}(:)))));
    legend(arrayfun(@(i) sprintf('机 %d', i), 1:n, 'UniformOutput', false), ...
        'Location', 'best');
end
xlabel('时间 (s)');
end

% ======================================================================
function animateScene(sim, cfg)
% 三维动画：n 架四旋翼（机体 + 四臂 + 旋翼）+ n 条绳索 + 长方体刚体负载。
nSteps = numel(sim.time);
n = cfg.vehicle.count;

% 显示坐标：把 z 取反（纵轴向上为正）
loadDisplay = sim.loadPositionLog;
loadDisplay(3, :) = -loadDisplay(3, :);
vehDisplay = squeeze(sim.vehiclePositionLog);
if n == 1
    vehDisplay = reshape(vehDisplay, 3, 1, nSteps);
end
vehDisplay(3, :, :) = -vehDisplay(3, :, :);

figure('Name', '多机协同吊运三维动画', 'Color', 'w', ...
    'Position', [100, 60, 1000, 720]);
ax = axes;
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

% ---- 坐标范围：按"起点 + 目标 + 最小跨度"确定 ----
% 不要用 min/max(整条轨迹) 自动适应：一旦超调较大，坐标轴会被拉到几十米，
% 0.1 m 量级的机体和负载会缩成一个点，完全看不出仿的是什么。
refPoints = [squeeze(vehDisplay(:, :, 1)), loadDisplay(:, 1), ...
    [cfg.target.position(1); cfg.target.position(2); -cfg.target.position(3)]];
limitLo = min(refPoints, [], 2) - cfg.visualization.axisPadding;
limitHi = max(refPoints, [], 2) + cfg.visualization.axisPadding;
span = max(limitHi - limitLo);
% ★ 坐标轴跨度必须与**负载实际尺寸和绳长**挂钩。
%   曾经默认 axisSpanMin = 2.00 m，导致 0.20 m 见方的负载只占画面的 10 %、
%   20 mm 的厚度只占 1 %，负载被压缩成一片看不出外形的薄板，用户完全看不出
%   仿的是什么。实测 0.55 m 时负载边长占 36 %、绳索占 64 %，才真正"看得见"。
payloadDiag = norm(cfg.payload.size);
spanMinEff = max([cfg.visualization.axisSpanMin, ...
    3.2 * payloadDiag, 1.4 * cfg.link.length]);
if span < spanMinEff
    mid = (limitLo + limitHi) / 2;
    limitLo = mid - spanMinEff / 2;
    limitHi = mid + spanMinEff / 2;
    span = spanMinEff;
end
xlim(ax, [limitLo(1), limitHi(1)]);
ylim(ax, [limitLo(2), limitHi(2)]);
zlim(ax, [limitLo(3), limitHi(3)]);
view(ax, 40, 20);
xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)'); zlabel(ax, '高度 (m)');

% ---- 几何模型（体系坐标）----
scaleV = cfg.visualization.vehicleScale;
armLength = cfg.vehicle.armLength * scaleV;
rotorRadius = cfg.vehicle.rotorRadius * scaleV;
[bodyVert, boxFaces] = boxVerticesFaces(cfg.vehicle.bodySize * scaleV);
[payloadVert, ~] = boxVerticesFaces( ...
    cfg.payload.size * cfg.visualization.payloadSizeScale);

% 把"负载是否够大"量化打印出来，便于判断显示参数是否合理。
% 经验判据：负载边长占跨度 > 20 % 才看得清外形，> 30 % 较舒适。
% ★ 这段必须放在 armLength 定义**之后**（以前放在前面，MATLAB 报
%   "函数或变量 'armLength' 无法识别"）。
%
% ★★★ 但**绝不能靠放大负载**来"看得清"：
%   挂点画在真实位置（由 payload.size 和 attachFractions 生成），四旋翼位置也由
%   veh = x0 + R0*rho_i - l*q_i 定死、绳索连到真实挂点。把负载按 k 倍画，
%   挂点就会落在"板面 1/(2k) 处"，**看起来像挂在板子中间** ——
%   用户正是据此判断"挂点不在边缘"。所以这里对 k≠1 明确报警。
if abs(cfg.visualization.payloadSizeScale - 1) > 1e-12
    warning('crazyflie_slung_visualization:PayloadScaleNotOne', ...
        ['payloadSizeScale = %.2f ≠ 1：负载被**放大显示**了，' ...
         '而挂点/绳索仍在真实位置 ⇒ 图上挂点会落在板面内部（约 %.0f%% 处），' ...
         '看起来像"挂在物品中间"。要看大请用图窗缩放，不要改这个系数。'], ...
        cfg.visualization.payloadSizeScale, ...
        100 / (2 * cfg.visualization.payloadSizeScale));
end
payloadSpanRatio = max(cfg.payload.size * cfg.visualization.payloadSizeScale) / span;
fprintf(['  [三维显示] 跨度 %.2f m | 负载边长占比 %.1f %% | ' ...
    '绳索占比 %.1f %% | 机臂占比 %.1f %%\n'], ...
    span, payloadSpanRatio * 100, cfg.link.length / span * 100, ...
    2 * armLength / span * 100);
if payloadSpanRatio < 0.20
    warning('crazyflie_slung_visualization:PayloadTooSmall', ...
        ['负载边长仅占坐标跨度的 %.1f %%，三维图里会看不清负载外形。' ...
         '请提高 cfg.visualization.payloadSizeScale（当前 %.2f）' ...
         '或减小 axisSpanMin（当前 %.2f）。'], ...
        payloadSpanRatio * 100, cfg.visualization.payloadSizeScale, ...
        cfg.visualization.axisSpanMin);
end

motorAngle = deg2rad([45, 135, 225, 315]);
motorBody = armLength * [cos(motorAngle); sin(motorAngle); zeros(1, 4)];

nRotorSeg = 24;
rotorUnit = [cos(linspace(0, 2 * pi, nRotorSeg + 1)); ...
             sin(linspace(0, 2 * pi, nRotorSeg + 1)); ...
             zeros(1, nRotorSeg + 1)];

% ---- 静态元素 ----
plot3(ax, cfg.target.position(1), cfg.target.position(2), ...
    -cfg.target.position(3), 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.85, 0.25, 0.10], 'Color', [0.85, 0.25, 0.10]);
plot3(ax, loadDisplay(1, :), loadDisplay(2, :), loadDisplay(3, :), '-', ...
    'Color', [0.75, 0.82, 0.92], 'LineWidth', 1.0);
drawTrajectories(ax, vehDisplay, n);

% ---- 动态元素：n 架四旋翼 + n 条绳索 + 负载 ----
linkLines = gobjects(n, 1);
bodyPatches = gobjects(n, 1);
armLines = gobjects(n, 4);
rotorLines = gobjects(n, 4);
bladeLines = gobjects(n, 4);
vehColors = lineColors(n);
for i = 1:n
    linkLines(i) = plot3(ax, nan, nan, nan, '-', 'Color', [0.30, 0.30, 0.30], ...
        'LineWidth', 2.0);
    bodyPatches(i) = patch(ax, 'Vertices', zeros(8, 3), 'Faces', boxFaces, ...
        'FaceColor', vehColors(i, :) * 0.75, 'FaceAlpha', 1.0, ...
        'EdgeColor', [0.05, 0.05, 0.05], 'LineWidth', 0.5);
    for j = 1:4
        armLines(i, j) = plot3(ax, nan, nan, nan, '-', ...
            'Color', [0.12, 0.12, 0.14], 'LineWidth', 2.6);
        rotorLines(i, j) = plot3(ax, nan, nan, nan, '-', ...
            'Color', [0.45, 0.60, 0.85], 'LineWidth', 0.9);
        if mod(j, 2) == 1
            bladeColor = [0.85, 0.20, 0.15];
        else
            bladeColor = [0.15, 0.35, 0.80];
        end
        bladeLines(i, j) = plot3(ax, nan, nan, nan, '-', ...
            'Color', bladeColor, 'LineWidth', 1.4);
    end
end
payloadPatch = patch(ax, 'Vertices', zeros(8, 3), 'Faces', boxFaces, ...
    'FaceColor', [0.55, 0.72, 0.92], 'FaceAlpha', 1.0, ...
    'EdgeColor', [0.15, 0.35, 0.65], 'LineWidth', 1.0);
% ★ 机体轴三色线：负载自转时若看不出朝向，视觉上会误以为负载变成了一个圆环。
%   （负载 yaw 误差较大时，30 s 内四角会扫出明显的圆形轨迹；
%   而负载厚度只有 0.02 m —— 侧视图上那圈就是扫掠轨迹。）
%   画上三条体轴后，自转方向与转速一眼可见，不会再被误读。
payloadAxisLines = gobjects(1, 3);
axisLen = 0.75 * max(cfg.payload.size);
payloadAxisColors = [0.85, 0.20, 0.15; 0.20, 0.60, 0.25; 0.20, 0.35, 0.85];
for j = 1:3
    payloadAxisLines(j) = plot3(ax, nan, nan, nan, '-', ...
        'Color', payloadAxisColors(j, :), 'LineWidth', 2.0);
end
attachMarkers = gobjects(n, 1);
for i = 1:n
    attachMarkers(i) = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.85, 0.25, 0.10], 'Color', 'none');
end

titleHandle = title(ax, sprintf('%d 机协同吊运仿真', n));

videoWriter = [];
if cfg.visualization.saveVideo
    videoWriter = VideoWriter(cfg.visualization.videoFile, 'MPEG-4');
    videoWriter.FrameRate = max(1, round(1 / (cfg.simulation.dt * ...
        cfg.visualization.animationStride)));
    open(videoWriter);
end

stride = cfg.visualization.animationStride;
for k = 1:stride:nSteps
    pl = loadDisplay(:, k);
    RlD = sim.loadRotationLog(:, :, k);
    RlD(3, :) = -RlD(3, :);                 % 体系 -> 显示系

    % 负载刚体
    set(payloadPatch, 'Vertices', (RlD * payloadVert + pl).');

    % 负载体轴（RlD 已是显示系），用于看清自转
    for j = 1:3
        tip = RlD * (axisLen * [j == 1; j == 2; j == 3]) + pl;
        set(payloadAxisLines(j), 'XData', [pl(1), tip(1)], ...
            'YData', [pl(2), tip(2)], 'ZData', [pl(3), tip(3)]);
    end

    % 挂点与绳索
    for i = 1:n
        pv = vehDisplay(:, i, k);
        RvD = sim.rotationLog(:, :, i, k);
        RvD(3, :) = -RvD(3, :);

        rhoInertial = sim.loadRotationLog(:, :, k) * cfg.payload.attachPoints(:, i);
        attach = pl + [rhoInertial(1); rhoInertial(2); -rhoInertial(3)];
        set(linkLines(i), 'XData', [pv(1), attach(1)], ...
            'YData', [pv(2), attach(2)], 'ZData', [pv(3), attach(3)]);
        if isfield(sim, 'takeoffModeLog')
            modeNow = sim.takeoffModeLog(k);
            if modeNow == 0 || modeNow == 1 || modeNow == 5
                set(linkLines(i), 'Color', [0.95, 0.55, 0.10], ...
                    'LineStyle', '--', 'LineWidth', 1.5);
            elseif modeNow == 2
                set(linkLines(i), 'Color', [0.95, 0.20, 0.10], ...
                    'LineStyle', '-', 'LineWidth', 2.2);
            else
                set(linkLines(i), 'Color', [0.30, 0.30, 0.30], ...
                    'LineStyle', '-', 'LineWidth', 2.0);
            end
        end
        set(attachMarkers(i), 'XData', attach(1), 'YData', attach(2), ...
            'ZData', attach(3));

        set(bodyPatches(i), 'Vertices', (RvD * bodyVert + pv).');

        spin = 2 * pi * cfg.vehicle.rotorSpinHz * sim.time(k);
        for j = 1:4
            motorDisplay = RvD * motorBody(:, j) + pv;
            set(armLines(i, j), 'XData', [pv(1), motorDisplay(1)], ...
                'YData', [pv(2), motorDisplay(2)], 'ZData', [pv(3), motorDisplay(3)]);

            ring = RvD * (rotorRadius * rotorUnit) + motorDisplay;
            set(rotorLines(i, j), 'XData', ring(1, :), 'YData', ring(2, :), ...
                'ZData', ring(3, :));

            th = spin + (j - 1) * pi / 2 + (i - 1) * pi / 4;
            blade = RvD * (rotorRadius * 0.92 * ...
                [-cos(th), 0, cos(th); -sin(th), 0, sin(th); 0, 0, 0]) ...
                + motorDisplay;
            set(bladeLines(i, j), 'XData', blade(1, :), 'YData', blade(2, :), ...
                'ZData', blade(3, :));
        end
    end

    % ★ 标题里显式给出负载偏航角，便于核对 yaw 参考与实际姿态是否一致。
    yawNow = atan2(sim.loadRotationLog(2, 1, k), sim.loadRotationLog(1, 1, k));
    modeText = '';
    if isfield(sim, 'takeoffModeLog')
        modeText = takeoffModeLabel(sim.takeoffModeLog(k));
    end
    slackNow = max(sim.ropeSlackLog(:, k));
    set(titleHandle, 'String', sprintf( ...
        ['t = %.2f s | %s | 负载高度 %.3f m | 负载偏航 %.1f deg | ' ...
         '绳余量 %.3f m | 张力 %.3f~%.3f N | 推力峰值 %.1f%%'], ...
        sim.time(k), modeText, -sim.loadPositionLog(3, k), rad2deg(yawNow), ...
        slackNow, min(sim.tensionLog(:, k)), max(sim.tensionLog(:, k)), ...
        max(sim.thrustPctLog(:, k))));
    drawnow;

    if ~isempty(videoWriter)
        writeVideo(videoWriter, getframe(gcf));
    end
end

if ~isempty(videoWriter)
    close(videoWriter);
end
end

% ======================================================================
function drawTrajectories(ax, vehDisplay, n)
colors = lineColors(n);
for i = 1:n
    plot3(ax, squeeze(vehDisplay(1, i, :)), squeeze(vehDisplay(2, i, :)), ...
        squeeze(vehDisplay(3, i, :)), '-', 'Color', colors(i, :) * 0.6 + 0.4, ...
        'LineWidth', 0.8);
end
end

function colors = lineColors(n)
% n 架四旋翼的区分色（生成式，不依赖工具箱）。
base = [0.85, 0.33, 0.10; 0.10, 0.55, 0.85; 0.20, 0.65, 0.30; ...
        0.60, 0.30, 0.75; 0.90, 0.65, 0.10; 0.15, 0.65, 0.65];
colors = zeros(n, 3);
for i = 1:n
    colors(i, :) = base(mod(i - 1, size(base, 1)) + 1, :);
end
end

function label = takeoffModeLabel(code)
switch code
    case 0
        label = 'SLACK 松弛';
    case 1
        label = 'TAKEUP 收紧';
    case 2
        label = 'TAUT_RAMP 软绷紧';
    case 3
        label = 'ACTIVE 协同运输';
    case 4
        label = 'LANDING_TAUT 受控下降';
    case 5
        label = 'LANDING_RELEASE 独立降落';
    otherwise
        label = 'UNKNOWN';
end
end

% ======================================================================
function y = unwrapLocal(x)
% 相位解卷绕（角度制）。等价于 Signal Processing Toolbox 的 unwrap(deg2rad(x))
% 再转回角度，但自写以保持本仿真"零工具箱依赖"。
% 判据：相邻样本跳变超过 180 度即认为跨过了 ±180 的割线，补偿 360 度。
x = x(:).';
y = x;
offset = 0;
for k = 2:numel(x)
    d = x(k) - x(k - 1);
    if d > 180
        offset = offset - 360;
    elseif d < -180
        offset = offset + 360;
    end
    y(k) = x(k) + offset;
end
end

% ======================================================================
function [verts, faces] = boxVerticesFaces(sizeVector)
% 长方体几何：8 个顶点（3x8，体系坐标）+ 6 个面的顶点索引。
s = sizeVector(:) / 2;
verts = [-s(1),  s(1),  s(1), -s(1), -s(1),  s(1),  s(1), -s(1); ...
         -s(2), -s(2),  s(2),  s(2), -s(2), -s(2),  s(2),  s(2); ...
         -s(3), -s(3), -s(3), -s(3),  s(3),  s(3),  s(3),  s(3)];
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end
