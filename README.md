# Geometric Control of Cable-Suspended Rigid Body

三架 Crazyflie 2.1 Brushless 协同吊运刚体负载的 MATLAB 仿真工程，基于 Lee 2014/2018 的几何控制框架。本分支使用八字轨迹、负载 yaw 参考和倾斜缆绳构型。

本版本重点处理：

1. 负载尺寸、挂点、惯量和相关控制参数保持一致；
2. 负载 yaw 使用低带宽闭环，不再关闭 yaw 反馈；
3. 缆绳允许倾斜，无人机不被强制放在负载正上方；
4. 小尺寸负载通过张力分配零空间的外张内部力降低无人机碰撞风险。

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

## 文件说明

| 文件 | 作用 |
|---|---|
| `crazyflie_slung_parameters.m` | 默认参数、负载尺寸派生量和合法性检查 |
| `crazyflie_slung_controller.m` | 负载位置/yaw 外环、张力分配、绳向控制和机体指令 |
| `crazyflie_slung_dynamics.m` | 负载、绳索和无人机动力学 |
| `crazyflie_slung_simulation.m` | 执行器模型、积分、日志和碰撞诊断 |
| `crazyflie_slung_reference.m` | 静态目标或八字轨迹参考 |
| `crazyflie_slung_demo.m` | 一键仿真和自检报告 |
| `crazyflie_slung_visualization.m` | 轨迹和三维构型可视化 |
| `crazyflie_slung_diagnose.m` | 数值异常和发散起点定位 |
| `ENGINEERING_LOG.md` | 公式、修改原因和调参记录 |

## 修改负载尺寸

`payload.size = [a; b; c]` 表示负载长、宽、高，单位为米。只修改尺寸时，程序会自动重建挂点、惯量、阻尼、负载姿态增益、悬停张力和初始推力：

```matlab
userCfg.payload.size = [0.12; 0.10; 0.06];  % [长; 宽; 高]
cfg = crazyflie_slung_parameters(userCfg);
sim = crazyflie_slung_simulation(cfg);
```

默认挂点模板为：

\[
\rho_1=(a/2,0,-c/2),\quad
\rho_2=(-a/2,b/2,-c/2),\quad
\rho_3=(-a/2,-b/2,-c/2).
\]

代码通过无量纲 `payload.attachFractions` 乘以 `payload.size` 生成实际挂点，因此尺寸变化后挂点仍位于负载上表面边界。显式提供 `userCfg.payload.attachPoints` 时，用户挂点优先，但必须为 `3 x n` 且列数等于无人机数量。

均质长方体惯量自动计算为：

\[
J_0=\operatorname{diag}\left(
\frac{m_0(b^2+c^2)}{12},
\frac{m_0(a^2+c^2)}{12},
\frac{m_0(a^2+b^2)}{12}
\right).
\]

同时更新 `payload.inertia`、`volume`、`surfaceArea`、`boundingRadius`、`rotationalDamping`、`loadController.kR` 和 `loadController.kOmega`。默认质量独立于尺寸；若要模拟同密度物品，省略 `mass` 并提供：

```matlab
userCfg.payload.size = [0.12; 0.10; 0.06];
userCfg.payload.density = 650;  % kg/m^3
cfg = crazyflie_slung_parameters(userCfg);
```

显式提供 `mass`、`inertia` 或 `rotationalDamping` 时，对应用户值优先。

## 负载 yaw 控制

负载 yaw 反馈默认开启：

```matlab
cfg.loadController.yawChannelEnabled = true;
```

八字轨迹默认 `figureEight.lockYaw = false`，因此期望 yaw 沿路径切线变化；设为 `true` 只会把参考 yaw 锁到 0，不会自动关闭反馈环。

控制器先计算负载 SO(3) 姿态误差和期望力矩 `Md`，再通过

\[
P=\begin{bmatrix}I&I&I\\\widehat\rho_1&\widehat\rho_2&\widehat\rho_3\end{bmatrix}
\]

把合力和合力矩分配到各根缆绳。yaw 力矩由缆绳的水平张力分量传递；默认目标带宽为 `0.80 Hz`，低带宽只表示避免激励摆动，不表示关闭 yaw。`yawChannelEnabled=false` 时才会在控制器中清零 yaw 力矩。

仿真输出包括 `loadYawLog`、`loadYawRefLog`、`loadYawErrorLog`、`steadyYawTrackingError` 和 `finalYawTrackingError`。无人机自身的 `headingSource` 只影响无人机绕推力轴的航向，不会关闭负载 yaw 环。

## 倾斜缆绳和碰撞避免

默认开启：

```matlab
cfg.link.allowTiltedCables = true;
```

控制器在最小范数张力解上加入内部力，并满足：

\[
P\mu_{\mathrm{internal}}=0.
\]

因此内部力不改变负载期望合力和合力矩，只改变各根缆绳的空间分布，使无人机从挂点正上方适度向外侧分开：

```matlab
cfg.allocation.outwardBiasFraction = 0.20;
cfg.allocation.outwardBiasMax = 0.12;
cfg.link.initialOutwardOffset = NaN;  % 自动设置
cfg.link.vehicleClearance = 0.02;
```

缆绳约束仍为：

\[
x_i=x_0+R_0\rho_i-l_iq_i,\qquad \|q_i\|=1.
\]

自检不再要求 `q_i=e_3` 或无人机位于负载正上方，而是检查绳长、张力正性、无人机间距、无人机与负载外接包络间隙和垂直净空诊断。

## 重要参数

| 参数 | 含义 |
|---|---|
| `payload.size` | 负载 `[长; 宽; 高]` |
| `payload.mass` / `payload.density` | 质量，或由密度和体积自动计算 |
| `payload.attachFractions` | 自动挂点模板 |
| `payload.attachPoints` | 负载坐标系中的实际挂点 `rho_i` |
| `payload.inertia` | 负载质心惯量矩阵 |
| `loadController.yawChannelEnabled` | 负载 yaw 反馈开关，默认 `true` |
| `link.allowTiltedCables` | 是否允许倾斜缆绳 |
| `allocation.outwardBiasFraction` | 外张内部力比例 |
| `allocation.outwardBiasMax` | 外张内部力上限 |
| `vehicle.collisionRadius` | 无人机碰撞包络半径 |
| `link.vehicleClearance` | 碰撞安全间隙 |

## 坐标和符号约定

- `q_i` 是从无人机指向负载的单位向量。
- `R_0` 和 `R_i` 是机体系到惯性系的旋转矩阵。
- `e_3=[0;0;1]` 指向重力方向，z 轴向下为正。
- 四旋翼实际作用力为 `-f_i R_i e_3`。
- 绳索只能受拉；张力非正表示绷紧缆绳模型不再适用。

## 修改后检查

1. 修改 `crazyflie_slung_parameters.m` 或通过 `userCfg` 覆盖参数。
2. 运行 `crazyflie_slung_demo('quick')`。
3. 检查 yaw 误差、最小张力、无人机间距和机体-负载间隙。
4. 再运行完整仿真并查看日志和图形。
5. 提交前运行 `git diff --check`。

本工程只包含本地代码和文档，不包含自动上传 GitHub 的脚本。
