# Geometric Control of Cable-Suspended Rigid Body

三架 Crazyflie 2.1 Brushless 协同吊运刚体负载的 MATLAB 仿真工程。模型和控制器对应 Lee 2014/2018 的几何控制框架，并保留 Crazyflie 推力执行器和角速度内环等效模型。

本版本重点处理五个问题：

1. 负载尺寸、挂点、惯量和相关控制参数保持一致（改 `payload.size` 后会自动重建派生量）；
2. 负载 yaw 使用**开启的低带宽闭环**，而不是关闭 yaw 反馈；
3. 缆绳允许倾斜，无人机不被强制放在负载正上方；
4. 小尺寸负载时通过张力分配零空间的外张内部力降低无人机碰撞风险，并按碰撞缺口自适应；
5. 加入"松弛绳起飞 → 收紧 → 软绷紧 → 协同运输 → 独立降落"的混合流程，
   并对交接瞬态做了逐项定量处理（见《地面起飞、绳索收紧与独立降落》）。

## 运行环境

- MATLAB R2021b 或更新版本
- 不依赖附加工具箱
- 建议运行前切换到本目录，不要把多个同名版本同时加入 MATLAB 路径。

## 快速运行

```matlab
% 先把本目录加入路径（换成你自己的路径）
addpath(pwd);
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
| `dump_sim_data.m` | （可选调试工具）跑一次仿真并把关键日志导出成 CSV + 诊断报告，便于用数据而非看图排查；不参与仿真，可直接删除 |
| `.gitignore` | 忽略 MATLAB 自动保存文件等（`*.asv` 等） |
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
（`omega_n = sqrt(kR/J0)` 直接随 `J0` 变）；更危险的是**阻尼比**——手填的 `kOmega`
与 `J0` 不匹配时，速率环等效极点 `kOmega/J0` 可能远高于环路采样频率
（`1/dt = 500 Hz`），姿态环会数值发散。

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

以当前默认参数（`payload.size = [0.08;0.06;0.05]`、`mass = 0.080 kg`）为例：

| 量 | 值 |
|---|---|
| `collisionRadius` | `armLength + rotorRadius` = 0.0695 m |
| `collisionRequiredSeparation` | `2x0.0695 + 0.05` = **0.1890 m** |
| `payloadNaturalSeparation` | 挂点 2 与挂点 3 距离 = **0.0600 m**（等于宽 b）|
| `outwardBiasScale` | **0.6825** |
| `outwardBiasFraction` | 0.20 缩小到 **0.1365** |
| 每机外张力（偏置 RMS）| 0.0906 N 缩小到 **0.0619 N** |
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

保持阻尼比不变，只把增益削到预算内。**默认尺寸不触发**
（`kR_x = 0.0578 < kRCap = 0.35·m0·g·rhoTyp/0.10 = 0.137`，其中 `rhoTyp = 0.05 m`）；
用户显式给出 `kR / kOmega` 时不削。派生量：`attitudeMomentCap`、`kRCap`、
`attitudeGainCapped`（是否被削）、`payload.horizontalArm`（即 `rhoTyp`，
取各挂点水平投影的最大值 = 0.05 m）。

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

`initialOutwardOffset` 只用于生成初始绳向；运行过程中缆绳方向由绳向动力学和绳向控制器决定。绷紧阶段（`TAUT_RAMP` / `ACTIVE`）缆绳长度约束保持：

$$
x_i=x_0+R_0\rho_i-l_iq_i,\qquad \|q_i\|=1.
$$

自检不再要求 `q_i = e_3` 或无人机水平投影必须在负载上方，而是检查绳长、张力正性、无人机间距、无人机与负载外接包络间隙，以及垂直净空诊断。

## 地面起飞、绳索收紧与独立降落

真实起飞不能从论文的绷紧约束直接开始。当前仿真默认
`cfg.takeoff.enabled = true`，状态机为：

```text
SLACK -> TAKEUP -> TAUT_RAMP -> ACTIVE
                                     |
                       LANDING_TAUT -> LANDING_RELEASE
```

> **命名提醒**：状态机里第三段的名字是 `ACTIVE`（见 `crazyflie_slung_simulation.m`
> 末尾 `modeCode` 的 `case`）。README 早期几版、`ENGINEERING_LOG.md` 和可视化图例
> 里写作 `TAUT_ACTIVE`，那只是叫法。**代码里 `strcmp(mode, ...)` 必须用 `'ACTIVE'`** ——
> 写错不会报错，只会让那个分支永远不命中（本项目为此踩过一次，见下）。

- `SLACK`：负载底面接触地面，三架无人机使用独立几何 PID 飞到安全高度；绳索长度由实际无人机位置计算，张力为 0。
- `TAKEUP`：无人机沿**悬停平衡绳向**（`link.takeupLinkUnitsBody`，由控制器求得）缓慢接近
  "挂点距离 = `link.length - preTensionSlack`"，用挂点距离和持续时间确认所有绳索接近绷直。
- `TAUT_RAMP`：在 `tensionRampTime` 内平滑建立协同控制输入（参考冻结在交接瞬间的实测状态），
  同时限制负载穿地。
- `ACTIVE`：进入 Lee 2014/2018 的绷紧缆绳模型和完整负载位置、姿态、yaw 控制。
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

### 交接参考剖面的设计（2026-09-28 修复，含实测依据）

交接段最容易搞错的是**参考的时序**。三条结论都是用转储数据（而非看曲线）得出的：

**① 张力建立期间负载**必然**贴地，参考不该动。**
本设计的悬停总张力恰好等于负载重量（`summary.hoverTensionByLink` 合计 = `m0*g`），
而 `tensionScale` 混合的 `uHover = m*g` 只抵无人机自重（**零张力**）。
所以 `tensionScale` 从 0 到 1 就是张力从 0 涨到 `m0*g`
⇒ **只有斜坡末端（`tensionScale ≈ 0.9~1.0`）负载才可能离地**。
实测：参考 1.5 s 内从 0.028 m 升到 0.350 m，而负载到 6.396 s 仍贴地
⇒ 位置误差在斜坡末端堆到 **0.3925 m**，之后负载才在 8.7 s 追上。
⇒ 正确顺序是**两段**：

```text
张力建立段 rampT（= tensionRampTime）  参考【冻结】在交接瞬间实测状态，误差恒 ≈ 0
抬升段     liftT（= referenceLiftTime）  参考用 5 次多项式剖面从冻结起点走到目标
```

**② 三个通道必须同源。**
`desired.position` 的导数要等于 `desired.velocity`、二阶导要等于
`desired.acceleration`，全部取自同一个 `smoothStep5WithDerivatives`。
只把位置做平滑、速度和加速度留 0，等于"参考自己在动、却要求速度为 0"，
位置环只能靠反馈硬追 ⇒ 负载明显滞后。
实测（旧写法，降落段同样问题）：参考 2.2 s 内从 0.350 降到 0.028 m，
负载只降到 0.146 m，**滞后 118 mm**，随后释放段不得不把负载"瞬移"到地面。

**③ 剖面的起点必须冻结。**
若锚点取当前实测值（`desired = 实测 + s·(目标 − 实测)`），误差只能按比例 `s` 释放；
负载不动时误差照样涨满 —— 这正是上一版的行为。所以交接瞬间要把
`tautStartPosition` / `tautStartVelocity` 快照下来（`landingStartPosition` 同理）。

参考剖面：交接后先冻结 `tensionRampTime`，再用 `referenceLiftTime` 走完到目标的位移。
峰值速度 ≈ `1.875·D/liftT`、峰值加速度 ≈ `5.77·D/liftT²`
（当前默认 `D ≈ 0.34 m`、`liftT = 3.0 s` ⇒ **0.214 m/s / 0.220 m/s²**；
旧写法把两者绑在同一个 1.5 s 斜披上时是 0.428 / 0.879）。
两端速度、加速度均为 0 ⇒ 与前面的地面段、后面的悬停段都 C² 连续。
用 `dump_sim_data.m` 导出的 `des_x/des_y/des_z` 可以逐点核对这条剖面。

**④ 顺带修掉的一个静默失效**（这类错误不报错，只让分支永不命中）：

```matlab
% ✗ 状态机里的名字是 'ACTIVE'，不是 'TAUT_ACTIVE' ⇒ 条件恒真
%   ⇒ 运输全程每一步都把 positionIntegral 清零，ki 完全失效
if ~(strcmp(mode, 'TAUT_ACTIVE') || strcmp(mode, 'LANDING_TAUT'))
    memory.positionIntegral = zeros(3, 1);
end
```

因为它是"清积分"的分支，失效后**没有任何报错或异常**，只是 `ki` 不再起作用。
修好之后 `ki` 才真正生效，所以巡航段行为会与之前略有不同，需要重新确认。
另外，交接抬升窗口内也必须继续清零积分器：那一段是**指令性瞬态**，
参考由前馈剖面给出，残余误差不是常值扰动；若让积分器累积，
抬升结束时它会带着"憋住的力"把负载顶过目标（新的过冲来源）。

### 收紧段绳向必须等于**悬停平衡绳向**（2026-09-28 修复）

交接瞬间的绳向误差 `eqi` 一开始就有大小之分，这不是"初始扰动"，而是**几何定义不一致**：

| | 收紧段（旧） | 悬停平衡（控制器期望） |
|---|---|---|
| 依据 | `takeoff.groundRadialOffset`（**一个统一值**） | 张力分配 + 零空间外张偏置 |
| 绳 1 | 16.4° | **10.93°** |
| 绳 2 | 16.4° | **15.27°** |
| 绳 3 | 16.4° | **15.27°** |

原因：悬停张力是 **2:1:1**（挂点质心偏离负载质心），受力大的那根绳更竖直
⇒ 控制器期望的平衡倾角是**三根各不相同**的。而旧写法用一个统一的水平偏移
把三机放到挂点外侧 ⇒ 三根绳倾角**完全相同**。
⇒ 交接瞬间 `eqi` 就有 5.6/1.2/1.0°，绳向环要在 0.2 s 内吞掉这个阶跃：
实测绳向误差冲到 **36°**、负载角速度 2.9 rad/s、机体速率与姿态指令**打到限幅**。

**修法**：收紧段的绳向改由**控制器自己的悬停平衡解**给出
（`simulation.m` 的 `equilibriumLinkUnitsBody()`：构造一个 `ex=0, ev=0, eR0=0` 的
`desired` 让控制器算一次分配，取回 `command.desiredLinkUnits`）：

```matlab
% 无人机 = 挂点 - (linkLength - preTensionSlack) * q_平衡   （q 由无人机指向负载）
cfg.link.takeupLinkUnitsBody = equilibriumLinkUnitsBody(state, cfg);
```

于是 **交接瞬间实际绳向 == 期望绳向 ⇒ `eqi ≡ 0`**（这是构造性成立，不是调参），
而 `q_id` 与 `tensionScale` 无关（`q_id = -mu_id/||mu_id||` 是尺度无关的），
所以整段 `TAUT_RAMP` 里 `eqi` 都保持 ≈ 0。

★ **顺带纠正一处概念混淆**：`groundRadialOffset` 是**避碰**量
（`boundingRadius + collisionRadius + vehicleClearance`），把它当**动力学**的绳向用
本来就不成立 —— 两个需求互不相关。地面段（`SLACK`）照旧用它，
收紧段改用平衡绳向。为免静默，`simulation.m` 会算一遍收紧段的机-负载最小水平距离：
低于**碰撞下限** `boundingRadius+collisionRadius` 才 `warning`；
低于"碰撞下限 + `vehicleClearance`"只打印一行提示（后者是**地面期**的额外余量，
机在负载上方 0.6 m 时不该再套用）。

### 地面必须给**水平摩擦**（2026-09-29 修复）

交接瞬态里还有一条持续扰动，来自一个纯粹的**建模缺口**：

```text
负载高度（转储实测）： 5.778 s → 0.0279 m   6.300 s → 0.0279 m
                       6.500 s → 0.0279 m   6.900 s → 0.0279 m
负载 x （转储实测）：  5.778 s → 0        6.900 s → −0.179 m
```

**负载一直贴在触地高度（0.0279 m = `groundZ − 半高`），却在水平方向滑了 0.18 m。**
原因是地面接触模型原来只有 z 向约束、**没有 xy 向摩擦** —— 等于"躺在无摩擦地面上"。
绳子的水平不平衡力（实测 ~0.07 N）把它拖走，挂点随之平移，
又反过来改变绳向、持续激励绳向环：绳向误差冲到 30°、负载角速度 2.7 rad/s、
速率/姿态指令打到限幅。

但真实情况下静摩擦上限 ≈ `μ·m0·g` = 0.24~0.39 N，是那个扰动的 **3~4 倍**
⇒ **负载根本不会滑**。所以原来的滑动是模型缺摩擦造出来的假象。

修法（`simulation.m` 地面接触块内）：接触期给库仑摩擦，水平减速度上限 `μ·g`，
并且 **μ 随法向力衰减** —— 绳把负载往上提 ⇒ 法向力 `N = m0·g − Σ T_i·q_i,z` 减小
⇒ 摩擦上限减小；张力涨到等于自重时 `N → 0`、摩擦自然消失。
这样离地瞬间是干净的，不需要额外的"释放判据"。

```matlab
'groundFrictionMu', 0.50, ...          % 负载触地时的库仑摩擦系数
'groundContactTolerance', 0.002, ...   % 判定"仍在地面接触"的竖直余量 [m]
```

### 独立控制器的积分器必须抗饱和（2026-09-29 修复）

同一轮还查出一个"看着小、但让判据变临界"的问题：收紧结束时绳长停在
**0.6200 m = `l − 0.0300`**，而 `epsilonOn = 0.030` —— **正好卡在绷紧判据的边界上**，
于是"绳是否绷紧"变成临界判断，TAKEUP 实际持续 3.28 s 而不是 `takeupDuration` 的 2.5 s。
（修好后实测 2.27 s，已回到 `takeupDuration` 附近。）

根因：`crazyflie_slung_independent_controller.m` 原来**无条件积分**，
而 SLACK/TAKEUP 是"从 0.3 m 外飞向目标"的大机动 ⇒ 积分器一路灌到限幅附近
（`Ki = 0.8`、限幅 0.20 ⇒ 最大 0.16 m/s² 的恒定力偏置）⇒ 停下来后无人机被稳稳压在
目标下方 **18 mm**，绳长也就到不了设计值 `l − preTensionSlack = 0.638`。

修法：加抗饱和 gate —— **位置误差大的时候不积分**：

```matlab
'independentIntegralGate', 0.050, ...   % 位置误差 < 它才积分 [m]
```

大机动段积分器不动（不影响跟踪），误差进入 gate 后才积分 ⇒ 稳态余差收干净，
收紧长度能落到设计值。交接处也顺手把该积分器清零（`TAUT_RAMP` 入口），
避免降落段重新启用独立控制器时带着起飞段的旧偏置。

### 交接点必须落在"绳刚好拉直"上（`preTensionSlack = 0`，2026-09-29 修复）

第三轮数据里最要命的一条：**模式切换那一拍，三架无人机同时瞬移 17 mm**（等价速率 8.5 m/s）。

```text
4.878 s（TAKEUP）      d = 0.6325 / 0.6326 / 0.6326   vehZ = −0.6702 / −0.6604 / −0.6684
4.880 s（TAUT_RAMP）   d = 0.6500 / 0.6500 / 0.6500   vehZ = −0.6873 / −0.6770 / −0.6853
                       ↑ 一拍（2 ms）内绳长被强制归零，无人机下移 17 mm
```

根因：**本模型只能表示绷紧的绳**（松弛段未实现）。而 `preTensionSlack = 0.012`
要求收紧末端保留 12 mm 余量 ⇒ 进入绷紧段时这 12 mm（加上独立控制器 ~5 mm 的稳态余差
= 17 mm）被一次性收掉 —— 数学上就是一次无人机位置瞬移。

⇒ 正确值就是 **0**：交接点是"绳刚好拉直、张力为零"，此时模型假设与实际一致、无跳变。
`parameters.m` 里对 `preTensionSlack` 的校验也从 `(0, L)` 放宽到 `[0, L)`；
`simulation.m` 在切换处会显式检查并**告警**（`takeupSnapWarn`，默认 10 mm），
不让这种隐性跳变悄悄过去。

> 残余量 ≈ 独立控制器的稳态余差（实测 ~7 mm，已低于 `takeupSnapWarn`）。要再压小它需要提高
> `takeoff.independentIntegralGain`（积分收敛时间 ≈ `2·Kp/Ki`，当前 0.8 ⇒ 约 10 s，
> 比 TAKEUP 的 2.4 s 长，所以收敛不完）。

### 地面还要约束**转动**（2026-09-29 补）—— 交接期"30° 绳向误差"的真源头

地面接触原来只约束 z（再加上后面补的水平摩擦），**完全不约束转动**。
于是交接期负载在地面上**打滚**：

```text
     t     负载 ωx      ωy      ωz    |ω|     负载高度
  4.850   −0.266   −1.502   +0.070   1.53    0.0302
  5.000   +0.571   +2.023   −0.033   2.10    0.0308
  5.050   +1.384   +2.286   −0.131   2.68    0.0279   ← 高度一直是触地值
```
躺在地面上的负载不可能以 ~3 rad/s 翻滚（那要求把一条边抬起来）。
而它不只是"不好看" —— 角速度经

```text
Md = -kR*eR0 - kOmega*eOmega0 + ...
```

直接进**张力分配**：`kOmega = [0.00276, 0.00403, 0.000339]`，`ω = 2.9 rad/s`
⇒ `|Md|` 可达 **0.012 N·m**；而力臂只有 ~0.04 m ⇒ 等效张力扰动
`0.012/(3*0.04) ≈ 0.10 N` —— **是绳 2 张力（0.20 N）的一半**
⇒ **期望绳向 `q_id` 被甩来甩去**，日志里 `eqi = cross(q_id, q_i)` 冲到 **30°**。

> ★★★ **但实际绳向根本没偏那么多。** 用无人机与负载位置**重建**实际绳向（`q_i = (x_0 + R_0 ρ_i − x_i)/‖·‖`，
> 只用到已经记录的 `sim.vehiclePositionLog` / `sim.loadPositionLog` / `sim.loadRotationLog`），
> 交接段"实际 q_i 与期望 q_id 的夹角"只有 **1~2°**（同一时刻日志里的 `ler` 是 29.6°）。
> ⇒ **那个 30° 主要是"被地面上的打滚甩出来的假象"**，不是绳真偏了。

物理上地面应该给恢复力矩：摩擦转矩上限 ≈ `μ*m0*g*rho_typ`，除以 `J0` 得
角减速度上限 `0.5*0.08*9.81*0.04 / 2.69e-4 ≈ 58 rad/s²`（滚转/俯仰），
所以触地期转动应当被**锁住**。修法（`simulation.m` 地面接触块内）：

```matlab
'groundRotationalBrake', 40.0, ...   % 触地时角减速度上限 [rad/s^2]
```

与水平摩擦一样**随法向力衰减**：张力把负载提起来时制动自然消失。

> ★ 诊断提醒：`sim.linkErrorLog` 是 `||cross(q_id, q_i)||`，**会被 `q_id` 自身的
> 抖动污染**。查"绳向到底偏没偏"时，除了看这个日志，还要用无人机与负载位置
> **重建实际绳向**再比一次 —— 两者不一致时，问题往往在"期望"那一侧。

### 交接瞬态：最终结果与残余（2026-09-29 收尾）

**排查过程留下的中间结论**（当时绳向峰值还在 ~30°，现已被下一条修复取代）：
表现为绳 2/3 从 16° 摆出到 36°/32°，而绳 1 几乎不动 —— 这是**低张力绳**效应：
绳 1 承担 2 倍张力（`|μ|` 0.399 vs 0.203 N），绳向模态的等效刚度 `~√(T/(m·l))`
也大约 2 倍 ⇒ 同样激励下绳 1 位移小、绳 2/3 位移大。
★ 但后来用位置重建发现，那 30° 里**大部分是"期望绳向 `q_id` 自己被甩"**，
真凶是下一条讲的"负载在地面上打滚"，见《地面还要约束转动》。

**修复后的最终指标**：

| 指标 | 修复前 | 现在 |
|---|---|---|
| 绳向误差峰值 | 29.6°（早期曾 36°）| **7.37°** |
| 负载角速度峰值 | 2.74 rad/s | **0.294 rad/s** |
| 交接期负载水平漂移 | −0.179 m | **−0.0004 m** |
| 模式切换位置瞬移 | 17.5 mm | 7.4 mm |
| 抬升段高度峰值（目标 0.350）| 0.3540（+4 mm）| **0.3506（+0.6 mm）** |
| 巡航段位置误差 | 0.11 mm | **0.3 mm** |
| 稳态绳向误差 | 0.096° | **0.081°** |

**残余（已量化，判定为可接受）**：

1. 峰值 7.37° 出现在交接后 0.066 s，主要是**绳 3 的初始偏向**。
   用位置重建实际绳向看（见上面的"诊断提醒"），
   绳 3 实际倾角 14.30° vs 平衡 15.27°、方位差 −3.31° ⇒ 真实夹角只有 **1.29°**
   ⇒ 日志里的 7° 仍有一部分来自 `q_id` 侧，但已是**单数量级**，不再是量级错误。
2. **张力建立段绳 1/3 有一个缓慢的"张开"**：倾角 12.0°→8.9°（绳1）、
   14.5°→22.8°（绳3），到 6.27 s 后收敛回平衡。`ler` 全程 ≤7° ⇒ 控制器跟得上，
   属**可接受的缓变**（对比修复前 36°/32° 的剧烈摆出）。
3. 模式切换残余瞬移 7.4 mm（原 17.5），来自独立控制器 ~7 mm 的稳态余差；
   已低于 `takeupSnapWarn = 10 mm`，不再告警。

**若还要继续压（未实施，按风险从低到高）**：
1. 提高 `linkController.kq` / `komega`（注意：本项目已知 3 机下 `kq` 有稳定上限，加大有风险）；
2. 缩短 `tensionRampTime`，减少处于"绳向环权限不足"区间的时间；
3. 让 `tensionScale` **只缩放与张力相关的部分**（`μ_i` 与绳向环项），
   `m·a_i` 这类惯性前馈不缩放 —— 机理上最对，但属控制器结构改动。

调整起飞流程时优先修改：

```matlab
cfg.takeoff.takeoffDuration    % 独立起飞时间
cfg.takeoff.takeupDuration     % 收紧时间
cfg.takeoff.tensionRampTime    % 软绷紧（张力建立）时间
cfg.takeoff.referenceLiftTime  % 张力建好后，参考抬升到目标的时间
cfg.takeoff.releaseSnapTolerance % 准许"落地点释放"的离地余量
cfg.takeoff.landingDuration    % 末段降落时间
cfg.takeoff.landingApproachFraction % 末段前一部分用绷紧动力学下降
cfg.takeoff.preTensionSlack    % 收紧末端保留余量
cfg.takeoff.groundRadialOffset % 地面阶段相对负载中心的安全外张距离
```

> `tensionRampTime` 只控制**张力建立**的快慢，不再同时决定参考走多远；
> 参考走多快由 `referenceLiftTime` 单独控制，两者解耦后互不干扰。

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
| `link.takeupLinkUnitsBody` | 收紧段绳向（负载体系，派生值；由 `simulation.m` 初始化时调用控制器求得）|
| `takeoff.preTensionSlack` | 收紧末端保留的绳长余量（默认 **0**，必须 0，见上文）|
| `takeoff.tensionRampTime` | 张力软建立时间（默认 1.5 s）|
| `takeoff.referenceLiftTime` | 张力建好后参考抬升到目标的时间（默认 3.0 s）|
| `takeoff.takeupSnapWarn` | 交接瞬间允许的绳长余量（= 位置瞬移量）告警阈值（默认 10 mm）|
| `takeoff.releaseSnapTolerance` | 准许"落地点释放"的离地余量（默认 5 mm）|
| `takeoff.groundFrictionMu` | 负载触地时的水平库仑摩擦系数（默认 0.50）|
| `takeoff.groundRotationalBrake` | 触地时角减速度上限（默认 40 rad/s²）|
| `takeoff.groundContactTolerance` | 判定"仍与地面接触"的竖直余量（默认 2 mm）|
| `takeoff.independentIntegralGate` | 独立控制器位置误差小于它才积分（抗饱和，默认 0.05 m）|
| `takeoff.*` | 其余地面接触、独立起飞/降落、绳索收紧和软张力过渡参数 |

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
6. **改过起飞/降落流程或地面接触后**，额外核对交接段这几个量：
   `summary.firstTautTime`、`sim.takeoffModeLog`（各段时长）、`sim.tensionScaleLog`、
   `sim.linkErrorLog`、`sim.loadBodyRateLog`、`sim.loadPositionLog`（看有没有水平漂移）。
   ★ 若要看得更细，直接跑 `dump_sim_data` 把关键日志导成 CSV + 诊断报告，
   用数据而不是看图来判断；`sim.linkErrorLog` 会被"期望绳向"自身的抖动污染，
   必要时用 `sim.vehiclePositionLog` / `sim.loadPositionLog` **重建实际绳向**再比一次。
7. 使用 `git diff --check` 检查格式。

## 免责声明

本工程是论文模型和 Crazyflie 执行器的数值仿真，不等同于真实飞行安全保证。修改负载尺寸、质量、缆绳长度或外张力后，必须重新检查张力正性、推力饱和、姿态误差和碰撞裕度。
