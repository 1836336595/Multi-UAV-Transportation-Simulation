# Geometric Control of Cable-Suspended Rigid Body

三架 Crazyflie 2.1 Brushless 协同吊运刚体负载的 MATLAB 仿真工程。模型和控制器对应 Lee 2014/2018 的几何控制框架，并保留 Crazyflie 推力执行器和角速度内环等效模型。

本版本重点处理四个问题：

1. 负载尺寸、挂点、惯量和相关控制参数保持一致；
2. 负载 yaw 使用开启的低带宽闭环，而不是关闭 yaw 反馈；
3. 缆绳允许倾斜，无人机不被强制放在负载正上方；
4. 小尺寸负载时通过张力分配零空间的外张内部力降低无人机碰撞风险。
5. 仿真加入“松弛绳起飞—收紧—软绷紧—协同运输—独立降落”的混合流程。

## 运行环境

- MATLAB R2021b 或更新版本
- 不依赖附加工具箱
- 建议运行前切换到本目录，不要把多个同名版本同时加入 MATLAB 路径。

## 快速运行

```matlab
cd('F:/workbuddy_doc/transporting/git')
report = crazyflie_slung_demo('quick');  % 快速自检，不绘图
report = crazyflie_slung_demo();         % 完整仿真、自检和可视化
```

建议先运行快速自检，再运行完整工况；仿真结果应以当前本地配置和 `summary` 输出为准。

## 文件说明

| 文件 | 作用 |
|---|---|
| `crazyflie_slung_parameters.m` | 默认参数、负载尺寸派生量和合法性检查 |
| `crazyflie_slung_controller.m` | 负载位置/yaw 外环、张力分配、绳向控制、推力和机体姿态指令 |
| `crazyflie_slung_independent_controller.m` | 松弛绳阶段的单机几何 PID 起飞/降落控制器（参考 v2） |
| `crazyflie_slung_dynamics.m` | 负载、绳索和无人机动力学 |
| `crazyflie_slung_simulation.m` | 执行器等效模型、积分、日志和碰撞诊断 |
| `crazyflie_slung_reference.m` | 静态目标或八字轨迹参考 |
| `crazyflie_slung_demo.m` | 一键仿真和自检报告 |
| `crazyflie_slung_visualization.m` | 轨迹、负载、无人机和缆绳可视化 |
| `crazyflie_slung_diagnose.m` | 发散起点和数值异常定位 |
| `ENGINEERING_LOG.md` | 公式、修改原因和调参记录 |

## 修改负载尺寸

只修改 `payload.size` 时，代码会自动重建与尺寸相关的派生参数：

当前默认值为 `payload.size = [0.08; 0.06; 0.05]` m；实际使用时以 `crazyflie_slung_parameters.m` 为准。

```matlab
userCfg = struct();
userCfg.payload.size = [0.12; 0.10; 0.06];  % [长; 宽; 高]，单位 m
cfg = crazyflie_slung_parameters(userCfg);
sim = crazyflie_slung_simulation(cfg);
```

`payload.size = [a; b; c]` 的含义是长、宽、高。默认挂点模板为：

$$
\rho_1=(a/2,0,-c/2),\quad
\rho_2=(-a/2,b/2,-c/2),\quad
\rho_3=(-a/2,-b/2,-c/2).
$$

实际代码通过无量纲 `payload.attachFractions` 乘以 `payload.size` 生成挂点，因此修改尺寸后挂点仍位于负载上表面边界。若显式提供 `userCfg.payload.attachPoints`，该自定义值优先，但必须为 `3 x n` 且列数等于无人机数量。

均质长方体惯量自动计算为：

$$
J_0=\operatorname{diag}\left(
\frac{m_0(b^2+c^2)}{12},
\frac{m_0(a^2+c^2)}{12},
\frac{m_0(a^2+b^2)}{12}
\right).
$$

同时自动更新 `payload.inertia`、体积、表面积、外接半径、转动阻尼、`loadController.kR`、`loadController.kOmega`、当前挂点几何对应的悬停张力和初始推力。

默认质量 `payload.mass` 是独立参数。如果希望材料密度不变、尺寸改变时质量随体积变化，可以省略 `mass` 并提供密度：

```matlab
userCfg.payload.size = [0.12; 0.10; 0.06];
userCfg.payload.density = 650;  % kg/m^3
cfg = crazyflie_slung_parameters(userCfg);
```

显式提供 `mass` 时，以 `mass` 为准；显式提供 `inertia` 或 `rotationalDamping` 时，对应参数也会保留用户值。

## Yaw 控制

负载 yaw 反馈默认开启：

```matlab
cfg.loadController.yawChannelEnabled = true;
```

控制器先计算负载 SO(3) 姿态误差和期望力矩 `Md`，再通过分配矩阵

$$
P=\begin{bmatrix}
I&I&I\\
\widehat\rho_1&\widehat\rho_2&\widehat\rho_3
\end{bmatrix}
$$

把合力和合力矩分配到各根缆绳。只要挂点不共线且缆绳有足够倾角，yaw 力矩就可以通过水平张力分量传递。yaw 使用较低带宽，是为了避免直接激励绳索摆动，不代表关闭 yaw 控制。

实际带宽（`parameters.m`）：

```matlab
wnLoad = 2*pi*[6.0; 6.0; 0.45];   % roll / pitch / yaw 的目标自然频率 [Hz]
zetaLoad = 0.90;
kRLoad     = diag(payloadInertia) .* wnLoad.^2;
kOmegaLoad = 2*zetaLoad*wnLoad .* diag(payloadInertia);
```

即 roll/pitch 6 Hz、**yaw 0.45 Hz**（约为前两者的 1/13）。低带宽的物理原因：
偏航力矩只能靠缆绳的**水平张力分量**传递，而水平张力同时会激励绳索摆动；
带宽越高、要求的偏航力矩越大、水平张力越大，越容易把绳摆激起来。
`designBandwidthHz` 与 `designDampingRatio` 是**尺寸派生量**的入口——
修改 `payload.size` 后 `kR`/`kOmega` 会按新的 `J0` 自动重算。

仿真记录：`sim.loadYawLog`、`sim.loadYawRefLog`、`sim.loadYawErrorLog`、`sim.summary.steadyYawTrackingError` 和 `sim.summary.finalYawTrackingError`。

负载 yaw 与无人机自身 yaw 是不同通道。`attitudeController.headingSource` 只决定无人机绕推力轴的机体航向参考，不会关闭负载 yaw 环。

### 负载姿态增益的三种设定方式

`cfg.loadController.kR / kOmega`（负载的 roll/pitch/yaw 三个轴各一个值）**默认是算出来的**，
不是手填的。来源是「目标带宽 + 阻尼比 + 负载惯量」：

```matlab
wnLoad   = 2*pi*[6.0; 6.0; 0.45];   % 目标带宽：roll/pitch 6 Hz，yaw 0.45 Hz
zetaLoad = 0.90;                    % 目标阻尼比
kR       = diag(J0) .* wnLoad.^2;   % 反解：omega_n = sqrt(kR/J0)
kOmega   = 2*zetaLoad*wnLoad .* diag(J0);
```

**为什么用这个写法**：`J0` 会被约掉 ⇒ 闭环带宽正好等于你指定的 `wnLoad`，
与负载多重多大无关。若把增益写死成常数，改尺寸后带宽会跟着漂
（历史教训：`kR=0.55` 配 `J0=2.69e-4` 时，`kOmega/J0` 达 1300 rad/s ≈ 207 Hz，
逼近 500 Hz 采样的奈奎斯特边界 ⇒ 姿态环数值发散）。

设定优先级（从高到低）：

| # | 方式 | 是否受张力预算上限削减 |
|---|---|---|
| 1 | `userCfg.loadController.kR / kOmega` | **不受** |
| 2 | `cfg.loadController.manualKR / manualKOmega` | **不受** |
| 3 | 自动（`designBandwidthHz` + `designDampingRatio` + 最终 `J0`）| 受 |

手动用法：

```matlab
cfg.loadController.manualKR     = [0.39; 0.39; 0.0043];
cfg.loadController.manualKOmega = [0.0189; 0.0189; 0.0027];
```

查 `cfg.loadController.attitudeGainManual` 可知当前是否用的是手动值。

## 倾斜缆绳与碰撞避免

默认配置：

```matlab
cfg.link.allowTiltedCables = true;
```

控制器在满足负载合力和合力矩的最小范数张力解上加入零空间内部力：

$$
P\mu_{\mathrm{internal}}=0.
$$

因此它不会改变负载的期望合力和合力矩，只会改变各根缆绳的空间分布，使无人机从挂点正上方适度向外侧分开。外张强度由以下参数控制：

```matlab
cfg.allocation.outwardBiasFraction = 0.20;
cfg.allocation.outwardBiasMax = 0.12;
cfg.link.initialOutwardOffset = NaN;  % 自动取机体包络 + 安全间隙
cfg.link.vehicleClearance = 0.05;
```

### 外张强度按碰撞缺口自适应

外张内部力和初始外张偏移**都不随负载尺寸变化**，但『无人机之间够不够开』这件事
**只在小负载时才需要外张**：自然机间距约等于挂点间距（随尺寸线性增长），
而要求机间距 `2*collisionRadius + vehicleClearance` 是常数。
因此 `parameters.m` 在尺寸派生阶段自动缩放二者：

```matlab
reqSep = 2*vehicle.collisionRadius + link.vehicleClearance;
natSep = payloadNaturalSeparation(payload.attachPoints);   % 挂点最小两两间距
allocation.outwardBiasScale = clamp((reqSep - natSep)/reqSep, 0, 1);
allocation.outwardBiasFraction = allocation.outwardBiasFraction * scale;
link.initialOutwardOffset     = link.initialOutwardOffset     * scale;
```

派生出三个可检查字段：`allocation.collisionRequiredSeparation`、
`allocation.payloadNaturalSeparation`、`allocation.outwardBiasScale`。
**用户若显式提供 `outwardBiasFraction` 或 `initialOutwardOffset`，则尊重用户值、不再缩放。**

以当前默认参数（`payload.size = [0.08;0.06;0.05]`、`mass = 0.033 kg`）为例：

| 量 | 值 |
|---|---|
| `collisionRadius` | `armLength + rotorRadius` = 0.0695 m |
| `collisionRequiredSeparation` | `2x0.0695 + 0.05` = **0.1890 m** |
| `payloadNaturalSeparation` | 挂点 2 与挂点 3 距离 = **0.0600 m**（等于宽 b）|
| `outwardBiasScale` | **0.6825** |
| `outwardBiasFraction` | 0.20 缩小到 **0.1365** |
| 每机外张力 | 0.0368 N 缩小到 **0.0251 N** |
| `initialOutwardOffset` | 0.1195 m 缩小到 **0.0816 m** |

而负载大到 `natSep >= reqSep`（本方参数集约 0.19 m）时 `scale` 归零，
外张机制**自动关闭**、张力裕度全部还给缆绳——避免大负载下『只付代价、不拿收益』。

> 该缩放是**线性启发式**（缺口 50% 就给 50% 的力），力与间距并非线性关系；
> 真正判据仍以 `summary.vehicleSeparationMargin`（实测最小中心距减要求值）为准。

### 负载放大后的自适应保护（2026-09-28 新增）

负载尺寸放大后，会出现几个**与尺寸强相关**的失稳源。`parameters.m` 在尺寸派生阶段
自动加了三道保护。它们的共同出发点是同一个事实：
**姿态力矩只能靠缆绳交付，而可用张力 `m0*g/n` 不随尺寸变。**

#### ① 姿态力矩的『张力预算上限』

`kR = J0 .* wn^2` 让期望力矩正比于 `J0`（尺寸的平方），而交付它靠各绳差分张力
`delta_mu ≈ Md/(n*rho)`。于是 `delta_mu/可用张力` **正比于尺寸**：

| 尺寸倍数 k | 1.0 | 2.0 | 2.5 | 3.0 | 4.0 |
|---|---|---|---|---|---|
| `delta_mu / 可用张力` | 14.7% | 29.5% | 36.8% | 44.2% | 58.9% |

⇒ 负载一大就把张力预算吃光 ⇒ 绳趋松弛、绳向环追不上 ⇒ 发散。

```matlab
momentCap = cfg.loadController.attitudeMomentBudget ...   % 默认 0.35
    * cfg.payload.mass * cfg.vehicle.gravity * rhoTyp;
cfg.loadController.kRCap = momentCap / cfg.loadController.attitudeMomentRefError;  % 0.10 rad
cfg.loadController.kR = min(cfg.loadController.kR, cfg.loadController.kRCap);
cfg.loadController.kOmega = 2 * zetaLoad .* sqrt(cfg.loadController.kR .* diag(cfg.payload.inertia));
```

保持阻尼比不变，只把增益削到预算内。**默认尺寸不触发**（`kR_x = 0.0238 < 0.0567`）；
用户显式给出 `kR / kOmega` 时不削。派生量：`attitudeMomentCap`、`kRCap`、
`attitudeGainCapped`（是否被削）、`payload.horizontalArm`（即 `rhoTyp`）。

#### ② 空中外张的『最小倾角地板』

`initialOutwardOffset` 原本只与机体参数有关（**尺寸的 0 次方**，等于 0.1195 m），
而力臂正比于尺寸 ⇒ 大负载时缆绳相对更竖直、水平张力可用量下降。加地板：

```matlab
cfg.link.initialOutwardOffset = max(cfg.link.initialOutwardOffset, ...
    cfg.link.minInFlightTiltRatio * rhoTyp);   % minInFlightTiltRatio 默认 0.20
```

#### ③ 增益尺寸防火墙

`kR / kOmega / kx / kv / ki` 必须都是 **3x1 列向量**，否则在参数阶段直接报错并
打印各量的实际尺寸。**教训**：曾把 `kOmega` 写成 `sqrt(...).'` 得到 1x3，
后续 `kOmega .* eOmega0` 因广播变成 3x3，错报在控制器第 125 行
`rhs6 = [R0.'*Fd; Md]`（vertcat 维度不一致）—— **报错点离病根很远，必须靠这道检查兜住**。

`initialOutwardOffset` 只用于生成初始绳向；运行过程中缆绳方向由绳向动力学和绳向控制器决定。`TAUT_RAMP` 和 `TAUT_ACTIVE` 阶段缆绳长度约束保持：

$$
x_i=x_0+R_0\rho_i-l_iq_i,\qquad \|q_i\|=1.
$$

自检不再要求 `q_i = e_3` 或无人机水平投影必须在负载上方，而是检查绳长、张力正性、无人机间距、无人机与负载外接包络间隙，以及垂直净空诊断。

## 地面起飞、绳索收紧与独立降落

真实起飞不能从论文的绷紧约束直接开始。当前仿真默认
`cfg.takeoff.enabled = true`，状态机为：

```text
SLACK -> TAKEUP -> TAUT_RAMP -> TAUT_ACTIVE
                                      |
                         LANDING_TAUT -> LANDING_RELEASE
```

- `SLACK`：负载底面接触地面，三架无人机使用独立几何 PID 飞到安全高度；绳索长度由实际无人机位置计算，张力为 0。
- `TAKEUP`：无人机缓慢接近 `link.length - preTensionSlack`，用挂点距离和持续时间确认所有绳索接近绷直。
- `TAUT_RAMP`：在 `tensionRampTime` 内平滑建立协同控制输入，同时限制负载穿地。
- `TAUT_ACTIVE`：进入 Lee 2014/2018 的绷紧缆绳模型和完整负载位置、姿态、yaw 控制。
- `LANDING_TAUT`：末段将负载平滑下降到地面；接触后进入 `LANDING_RELEASE`，无人机解除绳索约束并独立降落。

起飞阶段的坐标约定需要特别注意：动力学惯性系采用论文约定的
`e_3` 向下为正，而曲线图把 `-z` 显示为离地高度。因此地面为物理坐标
`z = groundZ`；机体中心最多到 `groundZ - vehicleGroundClearance`，负载中心的
接触高度按当前姿态下长方体的竖直包络实时计算。独立控制器中的位置/速度增益
先产生加速度，再乘以无人机质量得到力，避免把加速度增益误当成牛顿力导致起飞
振荡或穿地。

没有绳端拉力传感器时，仿真使用定位几何估计绳索状态：

```matlab
d_i = norm((x_0 + R_0*rho_i) - x_i);
q_i = ((x_0 + R_0*rho_i) - x_i) / d_i;
```

`epsilonOn`/`epsilonOff` 形成迟滞，`confirmTime` 防止噪声触发误切换。`sim.ropeDistanceLog`、`sim.ropeSlackLog`、`sim.takeoffModeLog` 和 `sim.tensionScaleLog` 可用于检查每次切换。绷紧阶段的 `sim.tensionLog` 是模型估计/指令张力，不是传感器实测值；真机应使用动捕/UWB 的无人机与负载位姿，并在有条件时增加绳端拉力传感器。

调整起飞流程时优先修改：

```matlab
cfg.takeoff.takeoffDuration    % 独立起飞时间
cfg.takeoff.takeupDuration     % 收紧时间
cfg.takeoff.tensionRampTime    % 软绷紧时间
cfg.takeoff.landingDuration    % 末段降落时间
cfg.takeoff.preTensionSlack    % 收紧末端保留余量
cfg.takeoff.groundRadialOffset % 地面阶段相对负载中心的安全外张距离
```

如果只想复现论文的“始终绷紧”模型，可设置 `cfg.takeoff.enabled = false`；此时无人机位置重新由 `x_i=x_0+R_0 rho_i-l_i q_i` 构造，不能用来验证地面起飞冲击。

## 重要参数

| 参数 | 含义 |
|---|---|
| `payload.size` | 负载 `[长; 宽; 高]` |
| `payload.mass` / `payload.density` | 质量，或由密度和体积自动计算 |
| `payload.attachFractions` | 自动挂点的无量纲模板 |
| `payload.attachPoints` | 负载坐标系中的实际挂点 `rho_i` |
| `payload.inertia` | 负载质心惯量矩阵 |
| `loadController.yawChannelEnabled` | 负载 yaw 反馈开关，默认 `true` |
| `link.allowTiltedCables` | 是否启用倾斜缆绳和外张内部力 |
| `allocation.outwardBiasFraction` | 外张内部力比例 |
| `allocation.outwardBiasMax` | 外张内部力上限 |
| `allocation.outwardBiasScale` | 外张缩放系数（按碰撞缺口自动算，0~1）|
| `allocation.collisionRequiredSeparation` | 要求的最小机间距（派生值）|
| `allocation.payloadNaturalSeparation` | 挂点最小两两间距（派生值）|
| `vehicle.collisionRadius` | 无人机碰撞包络半径 |
| `link.vehicleClearance` | 碰撞诊断安全间隙 |
| `loadController.attitudeMomentBudget` | 姿态力矩占张力预算的比例（默认 0.35）|
| `loadController.manualKR` | **手动指定负载姿态增益** `kR`（留空 `[]` = 自动）|
| `loadController.manualKOmega` | **手动指定负载角速度增益** `kOmega`（留空 `[]` = 自动）|
| `loadController.attitudeMomentRefError` | 折算力矩上限用的参考姿态误差（默认 0.10 rad）|
| `link.minInFlightTiltRatio` | 空中绳向最小倾角比例（默认 0.20）|
| `takeoff.*` | 地面接触、独立起飞/降落、绳索收紧和软张力过渡参数 |

## 代码约定

- `q_i` 定义为“从无人机指向负载”的单位向量。
- `R_0` 和 `R_i` 都是机体系到惯性系的旋转矩阵。
- 惯性系 `e_3=[0;0;1]` 指向重力方向，z 轴向下为正。
- 四旋翼实际作用力为 `-f_i R_i e_3`，推力方向由实际姿态决定。
- 绳索只能受拉；松弛阶段按 `tension = 0` 处理，绷紧阶段如果模型张力降到非正值，说明当前轨迹、尺寸或控制增益超出绷紧缆绳模型的适用范围。

## 修改和提交前检查

1. 修改 `crazyflie_slung_parameters.m` 或通过 `userCfg` 覆盖参数。
2. 检查尺寸、挂点列数和无人机数量一致。
3. 运行 `crazyflie_slung_demo('quick')` 做快速自检。
4. 检查 `summary` 中的 yaw 误差、最小张力、无人机间距和机体-负载间隙。
5. 再运行完整仿真并检查图形和日志。
6. 使用 `git diff --check` 检查格式；本工程不包含自动上传 GitHub 的脚本。

## 免责声明

本工程是论文模型和 Crazyflie 执行器的数值仿真，不等同于真实飞行安全保证。修改负载尺寸、质量、缆绳长度或外张力后，必须重新检查张力正性、推力饱和、姿态误差和碰撞裕度。
