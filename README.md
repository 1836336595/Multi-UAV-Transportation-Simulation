# Geometric Control of Cable-Suspended Rigid Body

三架 Crazyflie 2.1 Brushless 协同吊运刚体负载的 MATLAB 仿真工程。模型和控制器对应 Lee 2014/2018 的几何控制框架，并保留 Crazyflie 推力执行器和角速度内环等效模型。

本版本重点处理四个问题：

1. 负载尺寸、挂点、惯量和相关控制参数保持一致；
2. 负载 yaw 使用开启的低带宽闭环，而不是关闭 yaw 反馈；
3. 缆绳允许倾斜，无人机不被强制放在负载正上方；
4. 小尺寸负载时通过张力分配零空间的外张内部力降低无人机碰撞风险。

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

\[
\rho_1=(a/2,0,-c/2),\quad
\rho_2=(-a/2,b/2,-c/2),\quad
\rho_3=(-a/2,-b/2,-c/2).
\]

实际代码通过无量纲 `payload.attachFractions` 乘以 `payload.size` 生成挂点，因此修改尺寸后挂点仍位于负载上表面边界。若显式提供 `userCfg.payload.attachPoints`，该自定义值优先，但必须为 `3 x n` 且列数等于无人机数量。

均质长方体惯量自动计算为：

\[
J_0=\operatorname{diag}\left(
\frac{m_0(b^2+c^2)}{12},
\frac{m_0(a^2+c^2)}{12},
\frac{m_0(a^2+b^2)}{12}
\right).
\]

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

\[
P=\begin{bmatrix}
I&I&I\\
\widehat\rho_1&\widehat\rho_2&\widehat\rho_3
\end{bmatrix}
\]

把合力和合力矩分配到各根缆绳。只要挂点不共线且缆绳有足够倾角，yaw 力矩就可以通过水平张力分量传递。yaw 使用较低带宽，是为了避免直接激励绳索摆动，不代表关闭 yaw 控制。

仿真记录：`sim.loadYawLog`、`sim.loadYawRefLog`、`sim.loadYawErrorLog`、`sim.summary.steadyYawTrackingError` 和 `sim.summary.finalYawTrackingError`。

负载 yaw 与无人机自身 yaw 是不同通道。`attitudeController.headingSource` 只决定无人机绕推力轴的机体航向参考，不会关闭负载 yaw 环。

## 倾斜缆绳与碰撞避免

默认配置：

```matlab
cfg.link.allowTiltedCables = true;
```

控制器在满足负载合力和合力矩的最小范数张力解上加入零空间内部力：

\[
P\mu_{\mathrm{internal}}=0.
\]

因此它不会改变负载的期望合力和合力矩，只会改变各根缆绳的空间分布，使无人机从挂点正上方适度向外侧分开。外张强度由以下参数控制：

```matlab
cfg.allocation.outwardBiasFraction = 0.20;
cfg.allocation.outwardBiasMax = 0.12;
cfg.link.initialOutwardOffset = NaN;  % 自动取机体包络 + 安全间隙
cfg.link.vehicleClearance = 0.02;
```

`initialOutwardOffset` 只用于生成初始绳向；运行过程中缆绳方向由绳向动力学和绳向控制器决定。缆绳长度约束保持：

\[
x_i=x_0+R_0\rho_i-l_iq_i,\qquad \|q_i\|=1.
\]

自检不再要求 `q_i = e_3` 或无人机水平投影必须在负载上方，而是检查绳长、张力正性、无人机间距、无人机与负载外接包络间隙，以及垂直净空诊断。

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
| `vehicle.collisionRadius` | 无人机碰撞包络半径 |
| `link.vehicleClearance` | 碰撞诊断安全间隙 |

## 代码约定

- `q_i` 定义为“从无人机指向负载”的单位向量。
- `R_0` 和 `R_i` 都是机体系到惯性系的旋转矩阵。
- 惯性系 `e_3=[0;0;1]` 指向重力方向，z 轴向下为正。
- 四旋翼实际作用力为 `-f_i R_i e_3`，推力方向由实际姿态决定。
- 绳索只能受拉；如果张力降到非正值，说明当前轨迹、尺寸或控制增益超出绷紧缆绳模型的适用范围。

## 修改和提交前检查

1. 修改 `crazyflie_slung_parameters.m` 或通过 `userCfg` 覆盖参数。
2. 检查尺寸、挂点列数和无人机数量一致。
3. 运行 `crazyflie_slung_demo('quick')` 做快速自检。
4. 检查 `summary` 中的 yaw 误差、最小张力、无人机间距和机体-负载间隙。
5. 再运行完整仿真并检查图形和日志。
6. 使用 `git diff --check` 检查格式；本工程不包含自动上传 GitHub 的脚本。

## 免责声明

本工程是论文模型和 Crazyflie 执行器的数值仿真，不等同于真实飞行安全保证。修改负载尺寸、质量、缆绳长度或外张力后，必须重新检查张力正性、推力饱和、姿态误差和碰撞裕度。
