# 三机协同吊运仿真（Crazyflie 2.1 Brushless）

基于 **Lee 2014**（arXiv:1403.3684）与 **Lee 2018**（IEEE TCST）的几何控制框架，
用 MATLAB 实现**三架四旋翼协同吊运一块刚体负载**的完整仿真：起飞 → 沿 Gerono
双纽线（"八字"）绕行两个锥形障碍物 → 降落归位。

纯 MATLAB（R2021b）实现，**无附加工具箱依赖**。

---

## 目录结构

```text
.
├── README.md                    ← 本文件
├── ENGINEERING_LOG.md           ← 完整工程记录（约 160 KB，**改编参数前必读**）
├── crazyflie_slung_demo.m       ← 入口：跑仿真 + 打印配置与全部自检
├── crazyflie_slung_parameters.m ← 所有可调参数
├── crazyflie_slung_reference.m  ← 参考轨迹（起飞 / 环绕 / 降落）
├── crazyflie_slung_controller.m ← 几何控制器 (20)(21)(27)(36)-(40)
├── crazyflie_slung_dynamics.m   ← 刚体 + 绳索动力学 (5)-(8)
├── crazyflie_slung_simulation.m ← 主循环：执行器/速率环建模、积分、记录
├── crazyflie_slung_visualization.m ← 三维动画与结果图
├── crazyflie_slung_diagnose.m   ← 独立的发散定位脚本（不画图，只报告）
└── _verify_python/              ← 数值验证工具（见 §4）
```

> 历史版本不在本仓库内，统一放在 `../_archive/` 维护。

---

## 1. 快速开始

```matlab
>> crazyflie_slung_demo
```

跑完会打印**关键配置 + 完整自检结论**（每条都带实测值与判据），并弹出结果图。

### 切换工况

| 工况 | 怎么切 |
|---|---|
| **八字避障**（默认） | 保持 `cfg.referenceFcn = @(t) crazyflie_slung_reference(t, cfg)` |
| **仅定高悬停** | `cfg.referenceFcn = []`，并把 `cfg.obstacles.enabled` 置 `false` |

---

## 2. 已验证结果（八字避障，52 s）

```text
起飞段峰值 18.5 mm | 环绕段峰值 705.1 mm / 均值 43.3 mm | 降落段峰值 183.9 mm
降落终点残差 90.0 mm | 到锥最小净间隙 125.3 mm  PASS（判据 > 50 mm）
稳态绳索张力 [0.3928, 0.1962, 0.1962] N（= 挂点几何解出的 2:1:1）
绳索长度不变量 2.22e-16 m | 最小张力裕度 33.0% | 推力峰值 69.1%
最大机体角速度 3.949 rad/s | 状态有限 True | 坏步计数 0
```

---

## 3. 几个必须知道的结论

### 3.1 绳不是刚性连杆，但绷紧时力学等价
论文的动力学 (5)-(8)、分配 (13)(22)(23)、控制器 (27)(36)-(40) **一行都不用改**，
只多一条单边约束 `mu_i >= 0`。绳长不变量实测 **2.22e-16 m**（构造性成立）。

### 3.2 挂点几何决定静载分布（不是可选参数）
挂点取"**一边中点 + 对边两顶点**"时，挂点质心偏离负载质心 `a/3 = 0.0333 m`，
悬停张力因此是 **2:1:1**（`[0.3924, 0.1962, 0.1962] N`）而非三等分。
⇒ 机 1 推力需求提高、最小张力降低，**必须靠自检确认全程无松弛**。

### 3.3 `attitudeController.kR` 必须比直觉大得多（本构型 ×8）
不对称静载让各机受到的绳索反作用不再对称 ⇒ 机体姿态必须更"硬"，否则闭环发散。
实测 ×1 时静态悬停 1.202 s 发散；悬崖在 ×3~×5 之间，**×5~×20 是平台区**，取 ×8。

### 3.4 位置环带宽要取绳摆频率的 0.6~0.7 倍（反直觉）
`omega_p = sqrt(g/L) = 5.294 rad/s`。带宽越接近绳摆频率，位置环越会去追摆运动
产生的误差、把能量注入摆 ⇒ **KX=28（1.00×）反而发散最快（2.42 s）**。

### 3.5 ★★ 负载偏航通道**结构性不可控**
负载 yaw 是欠驱动自由轴（`kR0(3) = kOmega0(3) = 0`）。
**任何偏航反馈都会发散** —— 实测给 0.0183 的纯阻尼就冲到 1e11，
把负载放大到 1.0 m（论文量级）也一样。

**机理**：偏航力矩只能来自绳索的**水平分量**，而绳倾角又由负载自身状态决定
（`veh = x0 + R0·rho_i − l·q_i`）⇒ **"执行器"是"被控量"的函数** ⇒ 回路没有独立权限。

**能做的是**：
1. **参考端指向正确** —— 期望偏航 = 参考轨迹的**实际运动方向**
   （`psi = atan2(v_y, v_x)`，注意**不能用理想切向**，见 3.6）；
2. **实际端不动荡** —— 把**机体航向**与负载参考 yaw 解耦
   （`cfg.attitudeController.headingSource = 'worldX'`）。

### 3.6 ★★ 期望偏航 = 巡航段的**实际运动方向**（不是理想切向）

参考姿态第一轴取**参考轨迹的实际速度方向** `psi = atan2(v_y, v_x)`。

**为什么不能用"理想八字的解析切向"**：实际速度
`v = env'*p_hat + env*s*v_hat`，在包络过渡段（前/后 `blendTime` 秒）`|v|->0` 时
被 **`env'*p_hat`** 主导 ⇒ 真实运动是**径向**，而理想切向是**切向**。
实测两者在过渡段最大差 **93.1° / 162.7° / 178.3°**（τ=3.0 / 43.0 / 44.5 s）
—— 表现就是"偏航有时朝运动方向、有时正好相反"。

**为什么不能"直接取 atan2(v_y,v_x)"**：`|v|->0` 处方向会**翻转**
（实测 t=47.982 s 处 psi 由 +179.98° 一步跳到 0°，|v| 仅 1e-6）。

**最终实现**（`crazyflie_slung_reference.m`）：

```text
w   = v + epsV * t_hat_unit      % t_hat = 夹紧 tau 后的理想切向（恒不为零）
psi = atan2(w_y, w_x)            % epsV = 1e-4 m/s
```

* `|v| >> epsV` ⇒ `w ≈ v` ⇒ **就是实际运动方向**
* `|v| -> 0`  ⇒ `w ≈ epsV*t_hat` ⇒ 平滑接到过渡方向，**不跳 180°**

**实测：`|v| > 0.02` 时与真实运动方向的最大差 = 0.147°。**

**★ 起飞/降落段仍取 `psi = 0`**（不跟随）。原因有两条：
① 这两段近零速、方向病态；② 巡航两端 `|v|->0` 的极限方向**恰为 +x**，
与 `psi = 0` **天然衔接**，所以不会有跳变。
> ⚠ **曾把"跟随实际运动方向"推广到全阶段，结果出现安全性回归**：
> 起飞段水平位移方向是 149°，让参考 yaw 转 149° ⇒ 拧绳 ⇒
> **最小张力裕度 33.0% → 1.9%**（最小张力 0.0037 N，接近松弛）、推力峰值 69.1% → **100%**。
> **⇒ 已改回"仅巡航段跟随"。教训：参考更"真实"不等于闭环更好，推广前必须逐段检查安全裕度。**

**已知的固有代价**：巡航末端负载近乎静止时，速度方向本身无定义，
在约 0.02 s 窗口内会出现一次快速转向。这是"朝瞬时运动方向"这个定义的固有结果
（扫描 `epsV` 1e-4→2e-2 无干净解：即使 2e-2 跳变仍 14°/步、而方向偏差已涨到 77°）。
若要像素级连续，需对 `psi` 加时间低通（引入状态），当前不做。

**画图**：`crazyflie_slung_demo` 的**子图 6** 同时画**参考（虚线，`sim.loadYawRefLog`）
与 实际（实线）**，两条都 **unwrap 并对齐到同一支**（否则会出现"一条 +170、另一条 −190"
的假象），标题给出 rms 跟踪误差。
> ★ 参考 yaw 是在**仿真循环里记录**的（`sim.loadYawRefLog`），不是事后重调轨迹函数 ——
> 匿名句柄的多输出调用在部分 MATLAB 版本上会失败（这正是 `callReferenceFunction`
> 要 6→5→4→3 回退的原因），事后调用一旦失败会**静默退化**成"只画实际 yaw"。

### 3.7 ★ 机体航向锁定能显著减小负载偏航抖动
机体航向跟着负载参考 yaw 转 ⇒ 挂点方位转动 ⇒ 把绳"拧"起来 ⇒
反过来激励负载偏航（不可控、无法耗散）。
四旋翼的 **yaw 与推力解耦**，所以锁定机体航向的代价是**零**：

| `headingSource` | 负载偏航率抖幅(std) | 终端残差 |
|---|---|---|
| `'reference'`（论文原式） | 0.0498 rad/s | 94.0 mm |
| **`'worldX'`（默认）** | **0.0269 rad/s** | **90.0 mm** |

---

## 4. 数值验证（`_verify_python/`）

MATLAB 代码**同时有一份逐行对应的 Python 数值镜像**，所有定量结论都来自它；
MATLAB 侧另配一套静态检查。

| 工具 | 作用 |
|---|---|
| `_multi_mirror.py` | **主镜像**：逐行对应 MATLAB，数值结论的唯一来源 |
| `_lint_matlab.py` | MATLAB 静态检查 **7 类**：函数可见性 / 变量使用顺序 / `nargout` 分派 / `cfg` 路径 / 结构体字段 / 代码块配平 / **续行块语法** |
| `_test_lint_regress.py` | 上述检查的回归测试（保证"能抓到真错、不误报好码"） |
| `_test_reference_dispatch.py` | 参考函数多契约分派的审计 |
| `_test_obstacle_clearance.py` | 锥形障碍物间隙自检的回归测试（含"环绕性核对"） |

```bash
python _verify_python/_lint_matlab.py           # 静态检查
python _verify_python/_multi_mirror.py          # 跑数值镜像
python _verify_python/_test_lint_regress.py     # 回归测试
```

> ★ 为什么要镜像：开发环境**无法调用 MATLAB**，数值正确性只能靠一份逐行对应的
> Python 实现验证；而"文件私有函数 / `nargout` 语义 / 续行语法"这类
> **MATLAB 特有**的坑，镜像永远抓不到，只能靠 `_lint_matlab.py` 兜底。

---

## 5. 引用

```bibtex
@article{lee2018geometric,
  title={Geometric Control of Quadrotor UAVs Transporting a Cable-Suspended Rigid Body},
  author={Lee, Taeyoung},
  journal={IEEE Transactions on Control Systems Technology},
  year={2018}
}
@article{lee2014geometric,
  title={Geometric Control of Multiple Quadrotor UAVs Transporting a Cable-Suspended Rigid Body},
  author={Lee, Taeyoung},
  journal={arXiv:1403.3684},
  year={2014}
}
```

> 论文 PDF **未包含在本仓库中**（版权原因）。请从 IEEE Xplore / arXiv 获取。
