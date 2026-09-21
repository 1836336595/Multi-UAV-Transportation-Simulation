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
% ★ 轴框**不能**按负载轨迹自动适应。原因：负载从 z=0.45 m 爬到 0.35 m
%   （显示系即从 -0.45 到 -0.35），全程 z 极差 0.797 m，而负载边长只有
%   0.20 m、绳索 0.35 m —— 轴框被那段"爬升"撑到 0.8 m，负载只占 25 %。
%   改为以**四旋翼集群在稳态窗口的位置**为中心定框，负载爬升的那一段
%   会自然伸出框外（视觉上仍然完整可读，因为它是单调上升的一条线）。
%   子图只画前 trajectoryWindow 秒，避免 0.81 m 的初始瞬态占满画面。
subplot(2, 3, 1);
winEnd = min(nSteps, max(2, round(cfg.visualization.trajectoryWindow / ...
    cfg.simulation.dt) + 1));
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
title(sprintf('负载与 %d 架四旋翼轨迹（前 %.0f s，稳态局部放大）', n, ...
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
xlabel('时间 (s)'); ylabel('高度 (m)');
title('高度（负载 vs 各四旋翼）');
legend([{'负载'}, arrayfun(@(i) sprintf('机 %d', i), 1:n, 'UniformOutput', false), ...
    {'目标高度'}], 'Location', 'best');

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

% ---- (6) 负载偏航角（欠驱动自由轴）----
% ★ 这一格原来画"机 1 角速度指令 vs 实测"，但它与自检输出的
%   "最大机体角速度" 完全重复。负载偏航角才是本仿真最需要被看见的量：
%   它单调漂移、导致绳向误差缓增、也是三维图里"圆环"的来源。
subplot(2, 3, 6);
yawDeg = rad2deg(atan2(squeeze(sim.loadRotationLog(2, 1, :)), ...
    squeeze(sim.loadRotationLog(1, 1, :))));
% 用自写的 unwrapLocal 代替 MATLAB 的 unwrap()：后者属 Signal Processing
% Toolbox，本仿真声明不依赖任何工具箱（见文件头的说明）。
yawUnwrap = unwrapLocal(yawDeg);
plot(time, yawUnwrap, 'LineWidth', 1.6, 'Color', [0.55, 0.25, 0.65]);
hold on; grid on;
% 参考斜率：+0.057 deg/s（实测），用于说明"没有真稳态"
dtdy = [time(1), time(end)];
plot(dtdy, yawUnwrap(1) + 0.057 * (dtdy - time(1)), ':', ...
    'Color', [0.6, 0.6, 0.6], 'LineWidth', 1.0);
xlabel('时间 (s)'); ylabel('偏航角 (deg)');
title('负载偏航角（欠驱动自由轴，见 README §5.1）');
legend({'yaw', '参考斜率 0.057 deg/s'}, 'Location', 'best');
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
% ★ 八字工况下必须把**期望轨迹本身**也纳入坐标范围，
%   否则前 3 s 的起飞段（原地上升）会把轴框定得很小，
%   等负载飞到 ±1.5 m 时四旋翼直接跑出画面。
refPoints = [refPoints, referencePathPoints(cfg)];
% ★ 障碍物也要纳入，否则锥可能被裁掉一半。
if isfield(cfg, 'obstacles') && isfield(cfg.obstacles, 'enabled') ...
        && cfg.obstacles.enabled
    obsEdge = cfg.obstacles.positions + cfg.obstacles.radius;
    obsEdge = [obsEdge, cfg.obstacles.positions - cfg.obstacles.radius];
    refPoints = [refPoints, [obsEdge(1, :); obsEdge(2, :); ...
        -cfg.obstacles.baseZ * ones(1, size(obsEdge, 2))]];
end
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
%   挂点画在真实位置（x=±0.10 = 0.20 m 方板的边缘），四旋翼位置也由
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
% ★ 期望八字轨迹（细虚线）与两个锥形障碍物：对应论文 Fig. 3
%   "following a figure-eight curve around two obstacles represented by cones"。
%   画在动态元素之前，保证它们在最底层，不会挡住四旋翼和负载。
drawObstacles(ax, cfg);
drawDesiredPath(ax, cfg);
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
%   （本构型偏航是欠驱动自由轴，30 s 内自转约 60 度，四角扫出直径 0.28 m 的圆，
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

    % ★ 标题里显式给出负载偏航角。偏航是欠驱动自由轴（见 README §5.1），
    %   它会单调漂移；如果不显示出来，用户只会看到负载在转却不知道原因。
    yawNow = atan2(sim.loadRotationLog(2, 1, k), sim.loadRotationLog(1, 1, k));
    set(titleHandle, 'String', sprintf( ...
        ['t = %.2f s | 负载高度 %.3f m | 负载偏航 %.1f deg | ' ...
         '张力 %.3f~%.3f N | 推力峰值 %.1f%%'], ...
        sim.time(k), -sim.loadPositionLog(3, k), rad2deg(yawNow), ...
        min(sim.tensionLog(:, k)), max(sim.tensionLog(:, k)), ...
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

% ======================================================================
function pts = referencePathPoints(cfg)
% 采样期望八字轨迹（显示系，z 已取反），用于确定坐标轴范围。
pts = zeros(3, 0);
if ~isfield(cfg, 'referenceFcn') || isempty(cfg.referenceFcn) ...
        || ~isfield(cfg, 'figureEight')
    return;
end
nSamp = 400;
tGrid = linspace(0, cfg.simulation.duration, nSamp);
for k = 1:nSamp
    try
        [p, ~, ~, ~, ~, ~] = cfg.referenceFcn(tGrid(k));
        pts(:, k) = [p(1); p(2); -p(3)];
    catch
        return;   % 轨迹函数不可用时静默退回（不阻断可视化）
    end
end
end

% ======================================================================
function drawDesiredPath(ax, cfg)
% 画期望的八字轨迹：细虚线，环绕段用实色强调。
if ~isfield(cfg.visualization, 'showReferencePath') ...
        || ~cfg.visualization.showReferencePath
    return;
end
if ~isfield(cfg, 'referenceFcn') || isempty(cfg.referenceFcn)
    return;
end
nSamp = 600;
tGrid = linspace(0, cfg.simulation.duration, nSamp);
xyz = nan(3, nSamp);
for k = 1:nSamp
    try
        [p, ~, ~, ~, ~, ~] = cfg.referenceFcn(tGrid(k));
        xyz(:, k) = [p(1); p(2); -p(3)];
    catch
        xyz = nan(3, nSamp);
        break;
    end
end
if all(isnan(xyz(:)))
    return;
end
plot3(ax, xyz(1, :), xyz(2, :), xyz(3, :), '--', ...
    'Color', [0.55, 0.55, 0.60], 'LineWidth', 1.2);
% ★ 环绕段（即八字本体）用更醒目的绿色点线重复一遍，
%   让"绕八字"这件事在三维图里一眼可见。
t1 = cfg.figureEight.takeoffDuration;
t2 = t1 + cfg.figureEight.cruiseDuration;
inCruise = tGrid >= t1 & tGrid <= t2;
if any(inCruise)
    plot3(ax, xyz(1, inCruise), xyz(2, inCruise), xyz(3, inCruise), '-', ...
        'Color', [0.15, 0.62, 0.30], 'LineWidth', 1.8);
end
% 起点与环绕段入口标记
plot3(ax, xyz(1, 1), xyz(2, 1), xyz(3, 1), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.30, 0.60, 0.90], 'Color', [0.10, 0.30, 0.60]);
idxIn = find(inCruise, 1, 'first');
if ~isempty(idxIn)
    plot3(ax, xyz(1, idxIn), xyz(2, idxIn), xyz(3, idxIn), 's', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.95, 0.75, 0.15], 'Color', [0.55, 0.40, 0.05]);
end
end

% ======================================================================
function drawObstacles(ax, cfg)
% 两个锥形障碍物（论文 Fig. 3 "two obstacles represented by cones"）。
% 用 nSeg 边棱锥近似圆锥：顶点 + 底面圆环，纯显示，不参与动力学。
%
% ★ 面片必须显式构造，不能像折线那样"把点首尾拼起来"——
%   棱锥的顶点要参与每一个三角面，拼接顺序写错会画成一团乱麻。
if ~isfield(cfg, 'obstacles') || ~isfield(cfg.obstacles, 'enabled') ...
        || ~cfg.obstacles.enabled
    return;
end
nSeg = 24;
theta = linspace(0, 2 * pi, nSeg + 1);
% ★ 显示系是 z 向上，而物理系 z 向下为正，两者关系为 zDisplay = -zPhysical。
%   锥底在物理 z = baseZ（离地 hBase = -baseZ），锥尖在物理 z = baseZ - height
%   （即更高处）。换算到显示系：
%       底面高度 = -baseZ                （baseZ = 0 时即地面 0）
%       尖端高度 = -baseZ + height        （比底面**高** height）
%   ★ 旧代码写 zTop = zBase + height 且把尖端放在 zTop，
%     在 baseZ = 0 时会把锥**倒过来**画（尖端朝下扎进地里）。
zBase = -cfg.obstacles.baseZ;                 % 显示系里的锥底高度
zApex = zBase + cfg.obstacles.height;         % 显示系里的锥尖高度（更高）
nObs = size(cfg.obstacles.positions, 2);
for k = 1:nObs
    cx = cfg.obstacles.positions(1, k);
    cy = cfg.obstacles.positions(2, k);
    R = cfg.obstacles.radius;
    % 顶点表：1 = 锥尖，2..nSeg+2 = 底面圆环（最后一点与第 2 点重合）
    V = [cx, cy, zApex; ...
         cx + R * cos(theta(:)), cy + R * sin(theta(:)), zBase * ones(nSeg + 1, 1)];
    % 侧面三角面 + 底面扇形
    F = zeros(nSeg * 2, 3);
    for s = 1:nSeg
        F(s, :) = [1, s + 1, s + 2];
    end
    for s = 1:nSeg
        F(nSeg + s, :) = [nSeg + 2, s + 1, s + 2];
    end
    patch(ax, 'Vertices', V, 'Faces', F, ...
        'FaceColor', [0.92, 0.45, 0.20], 'FaceAlpha', 0.55, ...
        'EdgeColor', [0.65, 0.25, 0.08], 'LineWidth', 0.6);
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
