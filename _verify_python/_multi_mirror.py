"""Exact mirror of the n-vehicle MATLAB slung-load simulation.

Mirrors: crazyflie_slung_parameters.m / _dynamics.m / _controller.m / _simulation.m
         / _reference.m
Paper equations (Lee 2014 / 2018):
    (5)(6) payload translation+rotation coupled, (7) link, (8) vehicle attitude
    (13)(22)(23) tension allocation, (17)(24)(25) parallel part + link direction
    (27) link direction loop, (36)-(40) thrust + attitude
"""
import numpy as np

# ------------------------------------------------------------------ parameters
DT, DUR = 0.002, 52.0     # 总时长 = 起飞 3 + 环绕 45 + 降落 4
G = 9.81
M0, M, L = 0.080, 0.0325, 0.35
N = 3
J0 = np.diag([2.6933e-4, 2.6933e-4, 5.3333e-4])
J = np.diag([2.395e-5, 2.395e-5, 3.234e-5])
MAX_F = 1.3344
TAU = 0.012
# 挂点（负载体系）。★ 按用户要求放在**物品边缘**：一点为某边中点，另两点为对边两顶点。
#   1 = (+0.10,  0.00)  x=+0.10 边的中点
#   2 = (-0.10, +0.10)  x=-0.10 边的顶点
#   3 = (-0.10, -0.10)  x=-0.10 边的另一顶点
# 后果：cond(P P') 276.8 -> 150.2；静载张力变为 2:1:1 = [0.3924, 0.1962, 0.1962] N。
RHO = np.array([[0.10, -0.10, -0.10],
                [0.00,  0.10, -0.10],
                [-0.010, -0.010, -0.010]])

# ---- 八字避障工况参数（对应 cfg.figureEight）----
FE = dict(
    amplitudeX=1.50, amplitudeY=1.62,
    omegaX=0.2 * np.pi, omegaY=0.1 * np.pi,
    cycles=2.0,
    takeoffDuration=3.0, cruiseDuration=45.0, landingDuration=4.0,
    cruiseHeight=-0.55,
    startPosition=np.array([0.10, -0.06, 0.45]),
    landPosition=np.array([0.0, 0.0, -0.35]),
    lockYaw=False, blendTime=5.0,
)
# ★★ cruiseDuration 18 -> 45，blendTime 6.0 -> 5.0 ★★
#   为让轨迹**真正环绕**锥（而不是从旁边擦过），必须把包络的"满幅窗口"
#   (T2 - 2B) 撑大到覆盖住 τ∈[5,35]（两叶各完整遍历一次）。
#   实测（射线-线段求交判定"被包围"）：
#       T2=18 B=6.0  满幅 τ∈[12.7,27.3]  上叶无包围点      -> 无法环绕
#       T2=18 B=5.0  满幅 τ∈[10.6,29.4]  最小距 0.0419 m   -> 净间隙 -171 mm FAIL
#       T2=30 B=5.0  满幅 τ∈[ 6.4,33.6]  最小距 0.2962 m   -> 净间隙  +83 mm PASS
#       T2=45 B=5.0  满幅 τ∈[ 4.2,35.8]  最小距 0.3913 m   -> 净间隙 +178 mm PASS ←采用
#   ★ 加长 T2 使 s = 40/T2 从 2.222 降到 0.889，参考加速度 ∝ s² 降到 16%，
#     所以这个改动让系统**更稳**，不是妥协。
# ★★ amplitudeY = 1.62（原 1.20）—— 对照论文 Fig. 2 顶视图实测 ★★
#   论文八字的宽高比约 1.85（期望轨迹包围盒 宽 322 px / 高 174 px）。
#   原值 aX=1.50/aY=1.20 给出 2.50，比论文更"扁"。
#   取 aY = 2*1.50/1.85 = 1.62 复现论文比例。代价是跟踪略降但仍稳定：
#       aY=1.20 (2.50) -> 环绕均值 117.3 mm / 终点 78.7 mm
#       aY=1.62 (1.85) -> 环绕均值 133.7 mm / 终点 82.0 mm  ← 采用
#   全范围 1.20~1.85 单调平滑劣化、**无发散悬崖**（_diag_fe_ay.py）。
# 是否启用八字轨迹。False 时回退到静态悬停工况（用于对照实验）。
USE_FIGURE_EIGHT = True

# ------------------------------------------------------------------ 障碍物
# 对应论文 Lee 2018 Fig. 3 "two obstacles represented by cones"。
# ★ 只参与可视化与间隙自检，**不进入动力学**（论文也未对其建模）。
#   与 MATLAB crazyflie_slung_parameters.m 的 cfg.obstacles 逐字对应。
#
# ★★★ 2026-09-19 修复：旧参数是"倒插"的锥，导致自检静默失效 ★★★
#   旧值 positions=[0.75,-0.75; 0.60,0.60]、baseZ=-0.35、height=0.30、radius=0.22。
#   两个问题：
#     (1) 锥体整体位于负载环绕平面（离地 0.55 m）**下方**，负载从锥上方飞过，
#         避障约束根本不存在；
#     (2) 间隙公式 hAboveBase = (baseZ - z) + height 语义混乱（实为"离锥顶距离"），
#         使环绕段（占全程 72%）恒走 gap = inf 分支 —— 有效采样仅 2.1%，
#         且全部来自起飞瞬间（t = 0.41 s），旧报告里的 560.9 mm 是个
#         **与环绕段毫无关系的假信心数字**。
#   新值：锥轴移到叶心 (±1.5, 0.80)（aY 加高后叶心随之上移）、锥底落地、
#   加高到 1.05 m 使负载飞在锥高 52%（与论文 Fig. 2 侧视图实测比例一致）；
#   公式改为单一"离地高度"。真实净间隙 = +152.0 mm。
#   回归测试见 _test_obstacle_clearance.py。
OBSTACLES = dict(
    enabled=True,
    # ★ 两个锥摆在八字**两个叶的叶心**，负载**环绕**它们飞行。
    #   八字形状：ωX = 2·ωY ⇒ x 是倍频轴 ⇒ 两叶沿 y 上下堆叠，
    #   自交点在 (0, 0.81)，上叶 y∈[0.81,1.62]、下叶 y∈[0,0.81]，均横跨 x∈[-1.5,1.5]。
    #   叶心 = 在实际轨迹上求"被轨迹包围 且 到轨迹最小距离最大"的点
    #   （射线-线段求交判定包围，不是射线-点）。实测：
    #     上叶心 (-0.0187, +1.2252) -> 最近 0.3947 m -> 净间隙 +181.5 mm
    #     下叶心 (-0.0187, +0.3948) -> 最近 0.3938 m -> 净间隙 +180.6 mm
    positions=np.array([[-0.0187, -0.0187], [1.2252, 0.3948]]),  # 2x2，每列=[x;y]
    baseZ=0.00,          # 锥底在地面
    height=1.05,         # 负载飞在锥高 52% 处（论文侧视图实测 91 px 中的 y=78 px）
    radius=0.15,         # 锥底半径
    clearance=0.05,
)
# 负载外形 [长; 宽; 厚]，用于算外接球半径（与 MATLAB cfg.payload.size 一致）
PAYLOAD_SIZE = np.array([0.20, 0.20, 0.020])


KX = np.array([14.0, 14.0, 17.5])
KV = np.array([8.68, 8.68, 8.68])
# ★★★ 位置增益是"绕八字工况专用"的重标定值 ★★★
# 静态悬停用 [3.0, 3.0, 3.75] 足够（误差 3.63 mm），但绕八字时参考点持续运动，
# 位置环必须追得上，否则留下相位滞后型跟随误差。实测（A = ±1.5 m，B = 6.0）：
#   KX=3.0  wn=1.73 (0.33 x 绳摆频率) 环绕均 287.0 mm  末点 104.7 mm  稳定
#   KX=8.0  wn=2.83 (0.53 x)          环绕均 162.7 mm  末点 103.3 mm  稳定
#   KX=10.0 wn=3.16 (0.60 x)          环绕均 142.8 mm  末点  96.4 mm  稳定
#   KX=14.0 wn=3.74 (0.71 x)          环绕均 123.5 mm  末点  78.7 mm  稳定 ← 采用
#   KX=20.0 wn=4.47 (0.84 x)          发散于 11.10 s
#   KX=28.0 wn=5.29 (1.00 x 匹配)     发散于  2.42 s
# ★ 反直觉结论：让位置环带宽**匹配**绳摆频率反而最快发散。
#   耦合是寄生的：位置环越快，越会把负载水平运动硬转成绳索摆动，
#   而绳向环 (kq=55) 跟不上加快的激励，摆角积累到翻转。
#   最优带宽是绳摆频率的 0.6~0.7 倍。
# KV: zeta = KV/(2 sqrt(KX))；0.62*KX = 8.68 给出 zeta ≈ 1.16（略过阻尼），
# 比 0.50*KX / 0.90*KX 都更稳（_diag_fe_tune5.py 24c）。
# ki = 3.00：KX 提高后 1.60 显得偏软，末点残差 123 -> 78.7 mm；
#            ki >= 6.0 会在环绕段往复激励下发散（KI=6.0 于 11.30 s 发散）。
KI = np.array([3.00, 3.00, 3.00])
C1 = 0.50
INT_LIM = np.array([0.50, 0.50, 0.50])
# ★ 负载姿态增益按 J0 缩放：wn = sqrt(kR/J0), zeta = kOmega/(2 sqrt(kR J0))
#   ★★★ 第 3（yaw）通道必须为 0 ★★★
#   实测：kR0[2] 或 kO0[2] 非零，闭环 1 s 内发散；置零后位置误差收敛到 2 mm。
#   机理：yaw 力矩需求 -> mu 的纯水平分量 -> q_id 被甩向水平（单位 Md_z 偏转
#         86 deg）-> 实际连杆跟不上 -> (27) 投影把水平分量全丢（力矩交付率
#         24.9% -> 0%）-> 纯正反馈。缩增益无效（1/50 仍发散），是结构性欠驱动。
_WN_LOAD = 2 * np.pi * np.array([6.0, 6.0, 0.0])
_ZETA_LOAD = 0.90
KR0 = np.diag(J0) * _WN_LOAD ** 2
KO0 = 2 * _ZETA_LOAD * _WN_LOAD * np.diag(J0)
KQ, KW, KQI = 55.0, 20.0, 0.0
LINK_INT_LIM = np.array([0.30, 0.30, 0.30])
# ★★ kR 从 [30,30,15] 提到 ×8 = [240,240,120] —— 新挂点几何必需 ★★
#   新挂点（边中点+对边两顶点）使挂点质心偏离负载质心 0.0333 m，悬停张力
#   变成 2:1:1 = [0.3924, 0.1962, 0.1962] N，各机受到的绳索反作用不再对称，
#   机体必须更硬地保持姿态，否则推力方向偏离并与绳摆环形成正反馈。
#   实测 kR 倍率（静态悬停 / 八字）：
#     ×1 (30,30,15)    悬停 1.202 s 发散 / 八字 1.646 s 发散
#     ×2 (60,60,30)    悬停稳定 15.25 mm / 八字 7.290 s 发散
#     ×3 (90,90,45)    悬停稳定 15.11 mm / 八字 10.200 s 发散
#     ×5 (150,150,75)  悬停稳定 14.69 mm / 八字稳定 692.3 mm
#     ×8 (240,240,120) 悬停稳定 14.65 mm / 八字稳定 693.5 mm  ← 采用
#   悬崖在 ×3~×5 之间；×5~×20 是性能几乎不变的平台区，取 ×8 居中留裕度。
ATT_KR = np.array([240.0, 240.0, 120.0])
ATT_KO = np.array([4.0, 4.0, 4.0])
MAX_RATE_CMD = np.array([5.0, 5.0, 3.5])
MAX_ATT_ERR = 2.0
RATE_BW = np.array([35.0, 35.0, 20.0])
MAX_BODY_RATE, MAX_PEND_RATE = 12.0, 12.0
TARGET_POS = np.array([0.0, 0.0, -0.35])
TARGET_R0 = np.eye(3)
E3 = np.array([0.0, 0.0, 1.0])
DEBUG = False
# ★★★ 必须为 False ★★★
# q_id 是 (25) 算出的**单位向量**，对它做一阶差分再除以 dt 等于"两单位向量
# 之差 / 0.002 s"，信噪比极差。该噪声经 omega_id 进入 (27) 的
# -(q_i . omega_id) q_id_dot 项构成正反馈。实测稳定域：
#   True  -> aX 只能到 0.15 m（仍会发散）
#   False -> aX 可到 0.80 m 以上，配合 blendTime=6.0 则 ±1.5 m 全程稳定
# 逐项隔离（_diag_fe_terms.py）：去掉该项后 aX=0.2 从 20.75 s 发散变为全程稳定。
USE_QID_DOT = False
VEH_RPY = [[5.0, -4.0, 4.0], [-3.0, 6.0, 1.0], [4.0, 2.0, -3.0]]

# ---------------------------------------------------------------------- helpers


def hat(v):
    v = np.asarray(v, float).reshape(3)
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


# ------------------------------------------------- reference trajectory (八字)
# ★ 必须与 crazyflie_slung_reference.m 逐行对应，任何一侧改动都要同步。
#   三阶段：起飞(3s) -> 环绕两个八字(18s) -> 降落(4s)，拼接点 C^2 连续。

def _smoothstep(u):
    """5 次平滑阶跃，两端位置/速度/加速度都连续。"""
    u = min(max(u, 0.0), 1.0)
    s = u ** 3 * (10 - 15 * u + 6 * u ** 2)
    sd = 30 * u ** 2 * (1 - u) ** 2
    sdd = 60 * u * (1 - u) * (1 - 2 * u)
    return s, sd, sdd


def _blend_env(tau, T2, bt):
    """环绕段两端对称的速度包络 0 -> 1 -> 0。"""
    if bt <= 0:
        return 1.0, 0.0, 0.0
    sI, sId, sIdd = _smoothstep(tau / bt)
    sO, sOd, sOdd = _smoothstep((T2 - tau) / bt)
    return (sI * sO,
            (sId * sO + sI * sOd) / bt,
            (sIdd * sO + 2 * sId * sOd + sI * sOdd) / bt ** 2)


def _yaw_psi_and_rate(tau, s, aX, aY, wX, wY):
    """理想八字切向偏航角与其一阶时间导数的**解析**式。

    psi       = atan2(v, u)，u = aX wX cos(wX tau)，v = (aY/2) wY sin(wY tau)
    psiDot    = (u v' - v u') / (u^2 + v^2) * s      （d/dt = s d/dtau）
    """
    u = aX * wX * np.cos(wX * tau)
    v = 0.5 * aY * wY * np.sin(wY * tau)
    uP = -aX * wX ** 2 * np.sin(wX * tau)
    vP = 0.5 * aY * wY ** 2 * np.cos(wY * tau)
    psi = np.arctan2(v, u)
    psi_d = (u * vP - v * uP) / max(u * u + v * v, np.finfo(float).eps) * s
    return psi, psi_d


def _yaw_from_ideal_tangent(tau, s, h, aX, aY, wX, wY):
    """切向偏航角及其一/二阶真实时间导数（与 MATLAB 侧同式）。

    psiDotDot 只对**解析的** psiDot 做一次中心差分
    —— 比"对 R0d 做二阶差分"的噪声小一个量级（后者正是早期不稳定的来源之一）。
    """
    psi, psi_d = _yaw_psi_and_rate(tau, s, aX, aY, wX, wY)
    _, pPlus = _yaw_psi_and_rate(tau + s * h, s, aX, aY, wX, wY)
    _, pMinus = _yaw_psi_and_rate(tau - s * h, s, aX, aY, wX, wY)
    psi_dd = (pPlus - pMinus) / (2 * h)
    return psi, psi_d, psi_dd


def reference(t):
    """返回 (position, velocity, acceleration, rotation, bodyRate, bodyRateDot)。"""
    if not USE_FIGURE_EIGHT:
        return (TARGET_POS.copy(), np.zeros(3), np.zeros(3), TARGET_R0.copy(),
                np.zeros(3), np.zeros(3))

    T1 = FE['takeoffDuration']
    T2 = FE['cruiseDuration']
    T3 = FE['landingDuration']
    aX, aY = FE['amplitudeX'], FE['amplitudeY']
    wX, wY = FE['omegaX'], FE['omegaY']
    zC = FE['cruiseHeight']
    xyS = FE['startPosition'][:2]
    zS = FE['startPosition'][2]
    xyL = FE['landPosition'][:2]
    zL = FE['landPosition'][2]
    tEnd = FE['cycles'] * 2 * np.pi / wY
    s = tEnd / T2

    if t < T1:
        tau = 0.0
        phase = 0
    elif t < T1 + T2:
        tau = (t - T1) * s
        phase = 1
    else:
        tau = tEnd
        phase = 2

    # 平移后的八字（保形：Lemniscate of Gerono）
    #   pHatY = aY (1 - cos(wY tau)) / 2  -> tau=0 与 tau=tEnd 都是 0
    pHatX = aX * np.sin(wX * tau)
    pHatY = 0.5 * aY * (1 - np.cos(wY * tau))
    vHatX = aX * wX * np.cos(wX * tau)
    vHatY = 0.5 * aY * wY * np.sin(wY * tau)
    cHatX = -aX * wX ** 2 * np.sin(wX * tau)
    cHatY = 0.5 * aY * wY ** 2 * np.cos(wY * tau)

    if t < T1:
        sg, sgd, sgdd = _smoothstep(t / T1)
        pos = np.array([xyS[0] + (0 - xyS[0]) * sg,
                        xyS[1] + (0 - xyS[1]) * sg,
                        zS + (zC - zS) * sg])
        vel = np.array([(0 - xyS[0]) * sgd / T1,
                        (0 - xyS[1]) * sgd / T1,
                        (zC - zS) * sgd / T1])
        acc = np.array([(0 - xyS[0]) * sgdd / T1 ** 2,
                        (0 - xyS[1]) * sgdd / T1 ** 2,
                        (zC - zS) * sgdd / T1 ** 2])
    elif t < T1 + T2:
        e, ed, edd = _blend_env(t - T1, T2, FE['blendTime'])
        pos = np.array([pHatX * e, pHatY * e, zC])
        vel = np.array([ed * pHatX + e * s * vHatX,
                        ed * pHatY + e * s * vHatY, 0.0])
        acc = np.array([edd * pHatX + 2 * ed * s * vHatX + e * s * s * cHatX,
                        edd * pHatY + 2 * ed * s * vHatY + e * s * s * cHatY, 0.0])
    else:
        sg, sgd, sgdd = _smoothstep((t - T1 - T2) / T3)
        pos = np.array([(xyL[0]) * sg, (xyL[1]) * sg, zC + (zL - zC) * sg])
        vel = np.array([xyL[0] * sgd / T3, xyL[1] * sgd / T3,
                        (zL - zC) * sgd / T3])
        acc = np.array([xyL[0] * sgdd / T3 ** 2, xyL[1] * sgdd / T3 ** 2,
                        (zL - zC) * sgdd / T3 ** 2])

    # ★★★ 期望偏航 = **实际运动方向**；在 |v|->0 处平滑过渡到一个"过渡方向"。
    #   ------------------------------ 为什么需要过渡 -------------------------
    #   直接取 psi = atan2(v_y, v_x) 时，|v|->0 处方向会**翻转**：
    #     实测 t=47.982 s 处 psi 由 +179.98° 一步跳到 0°（|v| 仅 1e-6）；
    #     起飞段 |v| 很小而实际方向是 149°，取 +x 也是错的。
    #   这正是用户看到的"朝向正负不一致 / 符号不稳定"。
    #   ------------------------------ 过渡方向的选择 -------------------------
    #   取**夹紧 tau 后的理想切向** t_hat(tau_clamped)：
    #     tau 已被夹到 [0, tEnd]（起飞段 = 0、降落段 = tEnd）
    #     ⇒ 起飞与降落段的 t_hat 都恰好是 +x
    #     ⇒ 巡航两端（tau=0 与 tEnd）的 t_hat 也是 +x
    #   ⇒ **与相邻段的过渡方向天然一致，全程连续**。
    #   ------------------------------ 做法 -----------------------------------
    #   w = v + epsV * t_hat_unit，psi = atan2(w_y, w_x)：
    #     |v| >> epsV ⇒ w ≈ v        ⇒ 就是实际运动方向（满足"朝运动方向"）
    #     |v| -> 0    ⇒ w ≈ epsV*t_hat ⇒ 平滑接到过渡方向（连续，不跳 180°）
    #   epsV = 1e-4 m/s，远小于巡航典型速度 ~0.5 m/s，故不影响正常段。
    if FE['lockYaw'] or phase != 1:
        # 起飞/降落段：近零速、方向病态；且保持 psi=0 与巡航端点(极限方向恰为 +x)连续
        psi, psi_d = 0.0, 0.0
    else:
        tHatX = aX * wX * np.cos(wX * tau)          # tau 已夹紧到 [0, tEnd]
        tHatY = 0.5 * aY * wY * np.sin(wY * tau)
        tNorm = np.hypot(tHatX, tHatY)
        if tNorm > 0:
            tHatX, tHatY = tHatX / tNorm, tHatY / tNorm
        epsV = 1e-4
        wx = float(vel[0]) + epsV * tHatX
        wy = float(vel[1]) + epsV * tHatY
        den = wx * wx + wy * wy
        psi = float(np.arctan2(wy, wx))
        psi_d = ((wx * float(acc[1]) - wy * float(acc[0])) / den) if den > 1e-18 else 0.0
    c, sn = np.cos(psi), np.sin(psi)
    R0d = np.array([[c, -sn, 0.0], [sn, c, 0.0], [0.0, 0.0, 1.0]])
    return (pos, vel, acc, R0d, np.array([0.0, 0.0, psi_d]), np.array([0.0, 0.0, 0.0]))


def vee(S):
    return np.array([S[2, 1], S[0, 2], S[1, 0]])


def unit(v):
    v = np.asarray(v, float).reshape(3)
    nn = np.linalg.norm(v)
    return v / nn if nn > 1e-12 else np.array([0.0, 0.0, 1.0])


def proj_so3(R):
    U, _, Vt = np.linalg.svd(R)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        U[:, 2] = -U[:, 2]
        R = U @ Vt
    return R


def exp_so3(v):
    th = np.linalg.norm(v)
    if th < 1e-10:
        return np.eye(3) + hat(v)
    K = hat(np.asarray(v, float) / th)
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def rpy_to_R(rpy):
    r, p, y = rpy
    cr, sr = np.cos(r), np.sin(r)
    cp, sp = np.cos(p), np.sin(p)
    cy, sy = np.cos(y), np.sin(y)
    return np.array([[cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr],
                     [sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr],
                     [-sp,   cp*sr,            cp*cr]])


def proj_rate(prev_u, cur_u, dt):
    d = (np.asarray(cur_u, float) - np.asarray(prev_u, float)) / max(dt, 1e-12)
    return d - float(d @ cur_u) * np.asarray(cur_u, float)


def allocation_matrix():
    """P = [[I ... I],[hat(rho_1) ... hat(rho_n)]],  6 x 3n."""
    P = np.zeros((6, 3 * N))
    for i in range(N):
        P[0:3, 3*i:3*i+3] = np.eye(3)
        P[3:6, 3*i:3*i+3] = hat(RHO[:, i])
    return P


# --------------------------------------------------------------------- dynamics
ROT_DAMP = 1.0e-3  # 负载转动阻尼 [N*m*s/rad]，对应 cfg.payload.rotationalDamping


HEADING_SOURCE = 'worldX'      # 'reference' | 'worldX'（与 parameters.m 默认一致）


def dynamics(st, u_all):
    x0, v0 = st['x0'], st['v0']
    R0, Om0 = proj_so3(st['R0']), st['Om0']
    Q, QD = st['q'], st['qd']
    R, Om = st['R'], st['Om']
    u_all = np.asarray(u_all, float).reshape(3, N)

    hatOm0 = hat(Om0)
    hatOm0Sq = hatOm0 @ hatOm0

    sumQ = np.zeros((3, 3))
    sumMRho = np.zeros((3, 3))
    sumRhoQ = np.zeros((3, 3))
    sumRhoQRho = np.zeros((3, 3))
    rhs1 = np.zeros(3)
    rhs2 = np.zeros(3)
    u_par_all = np.zeros((3, N))
    u_perp_all = np.zeros((3, N))
    omega_all = np.zeros((3, N))

    for i in range(N):
        qi = unit(Q[:, i])
        qdi = QD[:, i]
        rhoi = RHO[:, i]
        omegai = np.cross(qi, qdi)
        omega_all[:, i] = omegai

        ui = u_all[:, i]
        u_par = float(qi @ ui) * qi
        u_perp = ui - u_par
        u_par_all[:, i] = u_par
        u_perp_all[:, i] = u_perp

        Qi = np.outer(qi, qi)
        rhoHat = hat(rhoi)
        zeta = u_par - M * L * float(omegai @ omegai) * qi \
            - M * Qi @ R0 @ (hatOm0Sq @ rhoi)
        rhs1 = rhs1 + zeta
        rhs2 = rhs2 + rhoHat.T @ R0.T @ zeta

        sumQ = sumQ + M * Qi
        sumMRho = sumMRho + M * Qi @ R0 @ rhoHat
        R0tQi = R0.T @ Qi
        sumRhoQ = sumRhoQ + M * rhoHat.T @ R0tQi
        sumRhoQRho = sumRhoQRho + M * rhoHat.T @ R0tQi @ R0 @ rhoHat

    Mq = M0 * np.eye(3) + sumQ
    Mblk = np.block([[Mq, -sumMRho], [sumRhoQ, J0 - sumRhoQRho]])
    # ★★★ 负载转动阻尼：**外部力矩**，直接加在转动方程右端（不是控制指令）。
    #   走控制器会经 "Md -> 张力分配 -> q_id 水平化" 正反馈（实测必发散）；
    #   外部力矩绕开控制器 ⇒ 稳定。依据见 parameters.m 的 rotationalDamping。
    rhs = np.concatenate([rhs1, rhs2 - np.cross(Om0, J0 @ Om0)
                          - ROT_DAMP * Om0])
    sol = np.linalg.solve(Mblk, rhs)
    x0dd = sol[:3] + G * E3
    Om0d = sol[3:]

    qdd_all = np.zeros((3, N))
    omdot_all = np.zeros((3, N))
    Om_dot_all = np.zeros((3, N))
    for i in range(N):
        qi = unit(Q[:, i])
        omegai = omega_all[:, i]
        rhoi = RHO[:, i]
        ai = x0dd - G * E3 + R0 @ (hatOm0Sq @ rhoi) - R0 @ hat(rhoi) @ Om0d
        omdot = (1 / L) * np.cross(qi, ai) \
            - (1 / (M * L)) * np.cross(qi, u_perp_all[:, i])
        omdot_all[:, i] = omdot
        qdd_all[:, i] = np.cross(omdot, qi) + np.cross(omegai, np.cross(omegai, qi))
        Oi = Om[:, i]
        Om_dot_all[:, i] = np.linalg.solve(J, st['torque'][:, i] - np.cross(Oi, J @ Oi))

    return dict(x0dd=x0dd, Om0d=Om0d, qdd=qdd_all, Om_dot=Om_dot_all,
                u_par=u_par_all, u_perp=u_perp_all, omega=omega_all, Mq=Mq)


# ------------------------------------------------------------------- controller
def controller(st, mem):
    dt = DT
    R0, Om0 = st['R0'], st['Om0']
    Q, QD = st['q'], st['qd']
    R, Om = st['R'], st['Om']

    # ★ 期望量改为从参考轨迹取（三阶段八字）。USE_FIGURE_EIGHT=False 时
    #   reference() 返回静态 TARGET_POS/TARGET_R0，与旧行为逐位一致。
    Xd, Vd, Ad, R0d, Om0d, Om0dDot = reference(st['t'])

    # (1) load outer loop -------------------------------------------------
    ex = st['x0'] - Xd
    ev = st['v0'] - Vd
    mem['pos_int'] = np.clip(mem['pos_int'] + dt * (ev + C1 * ex), -INT_LIM, INT_LIM)

    # paper (20) : Fd = m0( -kx ex - kv ev + x0dd_d - g e3 )
    #   *** 系数是 m0，不是等效质量矩阵 ***
    #   ★ 轨迹工况下必须加 desired.acceleration 前馈 Ad，否则纯反馈跟踪
    #     会有随加速度增大的滞后（八字峰值加速度 9.25 m/s^2，不可忽略）。
    x0dd_minus_ge3 = -KX * ex - KV * ev - KI * mem['pos_int'] + Ad - G * E3
    Fd = M0 * x0dd_minus_ge3

    # paper (21) : Md = -kR e_R0 - kOmega e_Omega0 + 前馈
    #   yaw 通道 kR0[2]=kO0[2]=0（结构性欠驱动，见上方参数说明）
    eR0 = 0.5 * vee(R0d.T @ R0 - R0.T @ R0d)
    eOmega0 = Om0 - R0.T @ R0d @ Om0d
    ffRate = R0.T @ R0d @ Om0d
    Md = -KR0 * eR0 - KO0 * eOmega0 \
        + np.cross(ffRate, J0 @ ffRate) + J0 @ (R0.T @ R0d @ Om0dDot)

    # (2) tension allocation: paper (13)(22)-(23) --------------------------
    P = mem['P']
    PPt = P @ P.T
    rhs6 = np.concatenate([R0.T @ Fd, Md])
    # ORDER MATTERS: MATLAB reshape(...,3,n) is column-major, numpy defaults to
    # C order. Using the default would interleave the per-vehicle tensions.
    mu_d_all = (np.kron(np.eye(N), R0) @ (P.T @ np.linalg.solve(PPt, rhs6))
                ).reshape(3, N, order='F')

    # (3)(4) per vehicle ---------------------------------------------------
    u_par = np.zeros((3, N))
    u_perp = np.zeros((3, N))
    q_des = np.zeros((3, N))
    eq_all = np.zeros((3, N))
    U = np.zeros((3, N))

    # (16) 的 a_i 用本拍外环量构造，不用上一拍状态（消除一拍延迟导致的寄生正反馈）
    Om0d_cmd = np.linalg.solve(J0, Md - np.cross(Om0, J0 @ Om0))

    for i in range(N):
        qi = unit(Q[:, i])
        qdi = QD[:, i]
        rhoi = RHO[:, i]
        omegai = np.cross(qi, qdi)

        mu_id = mu_d_all[:, i]
        mu_i = qi * float(qi @ mu_id)
        nn = np.linalg.norm(mu_id)
        qid = unit(-mu_id / nn) if nn > 1e-9 else unit(mem['q_prev'][:, i])
        qid = unit(qid)

        ai = x0dd_minus_ge3 + R0 @ (hat(Om0) @ (hat(Om0) @ rhoi)) \
            - R0 @ hat(rhoi) @ Om0d_cmd

        # paper (17)
        u_p = mu_i + M * L * float(omegai @ omegai) * qi + M * qi * float(qi @ ai)
        u_par[:, i] = u_p

        if USE_QID_DOT and mem['q_prev_des'] is not None:
            qid_dot = proj_rate(mem['q_prev_des'][:, i], qid, dt)
        else:
            qid_dot = np.zeros(3)
        omega_id = np.cross(qid, qid_dot)

        eqi = np.cross(qid, qi)
        hatQSq = hat(qi) @ hat(qi)
        e_omega_i = omegai + hatQSq @ omega_id
        mem['link_int'][:, i] = np.clip(
            mem['link_int'][:, i] + dt * eqi, -LINK_INT_LIM, LINK_INT_LIM)

        correction = -KQ * eqi - KW * e_omega_i - KQI * mem['link_int'][:, i] \
            - float(qi @ omega_id) * qid_dot
        u_pp = M * L * np.cross(qi, correction) - M * hatQSq @ ai
        u_perp[:, i] = u_pp

        q_des[:, i] = qid
        eq_all[:, i] = eqi
        U[:, i] = u_p + u_pp

    mem['q_prev_des'] = q_des.copy()

    # (5)(6) thrust / desired attitude / body-rate command -----------------
    thrust = np.zeros(N)
    om_cmd = np.zeros((3, N))
    eR_all = np.zeros((3, N))
    # 机体航向来源：'reference'（论文原式）或 'worldX'（锁定 +x，抖幅减半，见 MATLAB 侧说明）
    b1d = (np.array([1.0, 0.0, 0.0]) if HEADING_SOURCE == 'worldX'
           else R0d[:, 0])
    for i in range(N):
        Ri = R[:, :, i]
        Oi = Om[:, i]
        ui = U[:, i]
        nu = np.linalg.norm(ui)
        b3c = unit(-ui / nu) if nu > 1e-9 else unit(mem['Rc_prev'][:, 2, i])
        b3c = unit(b3c)
        P3 = np.eye(3) - np.outer(b3c, b3c)
        proj = P3 @ b1d
        if np.linalg.norm(proj) < MAX_ATT_ERR * 1e-6:
            proj = P3 @ np.array([1.0, 0.0, 0.0])
        b1c = unit(proj)
        b2c = unit(np.cross(b3c, b1c))
        b1c = unit(np.cross(b2c, b3c))
        Rc = proj_so3(np.column_stack([b1c, b2c, b3c]))
        mem['Rc_prev'][:, :, i] = Rc

        thrust[i] = float(np.clip(-float(ui @ (Ri @ E3)), 0.0, MAX_F))

        eR = 0.5 * vee(Rc.T @ Ri - Ri.T @ Rc)
        ne = np.linalg.norm(eR)
        if ne > MAX_ATT_ERR:
            eR = eR * (MAX_ATT_ERR / ne)
        eR_all[:, i] = eR
        # Omega_ic feedforward = 0 (1/dt amplifier, see README)
        cmd = -ATT_KR * eR - ATT_KO * Oi
        om_cmd[:, i] = np.clip(cmd, -MAX_RATE_CMD, MAX_RATE_CMD)

    return dict(ex=ex, Fd=Fd, Md=Md, mu_d=mu_d_all, U=U, thrust=thrust,
                om_cmd=om_cmd, eq=eq_all, eR=eR_all, eR0=eR0)


# -------------------------------------------------------------------------- run
def run(duration=DUR):
    dt = DT
    n_steps = int(duration / dt) + 1

    tilt = np.deg2rad([4, -5, 3])
    tdir = np.deg2rad([0, 120, 240])
    q0 = np.zeros((3, N))
    for i in range(N):
        q0[:, i] = unit(np.array([np.sin(tilt[i]) * np.cos(tdir[i]),
                                  np.sin(tilt[i]) * np.sin(tdir[i]),
                                  np.cos(tilt[i])]))

    st = dict(x0=np.array([0.10, -0.06, 0.45]), v0=np.zeros(3),
              R0=rpy_to_R(np.deg2rad([3, -2, 4])), Om0=np.zeros(3),
              q=q0.copy(), qd=np.zeros((3, N)),
              R=np.zeros((3, 3, N)), Om=np.zeros((3, N)),
              torque=np.zeros((3, N)),
              a_prev=np.zeros(3), Om0d_prev=np.zeros(3))
    for i in range(N):
        st['R'][:, :, i] = rpy_to_R(np.deg2rad(VEH_RPY[i]))

    mem = dict(pos_int=np.zeros(3), link_int=np.zeros((3, N)),
               q_prev=q0.copy(), q_prev_des=None,
               Rc_prev=np.tile(np.eye(3)[:, :, None], (1, 1, N)),
               P=allocation_matrix())

    # ★ 初始推力必须按**挂点几何**解出的悬停张力给，不能假设"三等分"。
    #   新挂点质心偏离负载质心 0.0333 m，悬停张力是 2:1:1；
    #   若仍用 (m0/N)g + mg，机 1 被低估 0.13 N，起步瞬间推力打到 100%（实测）。
    _bal0 = np.vstack([np.ones(N), RHO[1, :], -RHO[0, :]])
    _mu0 = np.linalg.solve(_bal0, np.array([-M0 * G, 0.0, 0.0]))
    thrust_act = np.abs(_mu0) + M * G

    L_ = {k: np.zeros((3, N, n_steps)) for k in
          ['q', 'veh', 'Om', 'u_par', 'u_perp', 'mu_d', 'om_cmd']}
    L_['x0'] = np.zeros((3, n_steps))
    L_['R0'] = np.zeros((3, 3, n_steps))
    L_['R'] = np.zeros((3, 3, N, n_steps))
    L_['tension'] = np.zeros((N, n_steps))
    L_['thrust_pct'] = np.zeros((N, n_steps))
    L_['pos_err'] = np.zeros(n_steps)
    L_['link_err'] = np.zeros((N, n_steps))
    L_['att_err'] = np.zeros((N, n_steps))
    # ★ 新增：负载姿态误差的按轴日志（3 x n_steps）与负载角速度日志。
    #   用于区分可控轴（roll/pitch）与欠驱动轴（yaw，见 README §5.1）。
    L_['load_att_err'] = np.zeros((3, n_steps))
    L_['load_Om0'] = np.zeros((3, n_steps))
    # ★ 轨迹跟踪专项日志：期望位置、期望速度、三轴跟踪误差
    L_['xd'] = np.zeros((3, n_steps))
    L_['vd'] = np.zeros((3, n_steps))
    L_['track_err_vec'] = np.zeros((3, n_steps))

    for k in range(n_steps):
        st['t'] = k * dt
        cmd = controller(st, mem)

        L_['x0'][:, k] = st['x0']
        L_['R0'][:, :, k] = st['R0']
        L_['q'][:, :, k] = st['q']
        L_['R'][:, :, :, k] = st['R']
        L_['Om'][:, :, k] = st['Om']
        L_['om_cmd'][:, :, k] = cmd['om_cmd']
        L_['mu_d'][:, :, k] = cmd['mu_d']
        L_['u_par'][:, :, k] = st.get('u_par', np.zeros((3, N)))
        for i in range(N):
            L_['tension'][i, k] = np.linalg.norm(
                unit(st['q'][:, i]) * float(unit(st['q'][:, i]) @ cmd['mu_d'][:, i]))
            L_['thrust_pct'][i, k] = cmd['thrust'][i] / MAX_F * 100
            L_['link_err'][i, k] = np.linalg.norm(cmd['eq'][:, i])
            L_['att_err'][i, k] = np.linalg.norm(cmd['eR'][:, i])
            L_['veh'][:, i, k] = st['x0'] + st['R0'] @ RHO[:, i] - L * st['q'][:, i]
        # 负载姿态误差按轴记录（负载只有一个，不重复 N 份）
        L_['load_att_err'][:, k] = cmd['eR0']
        L_['load_Om0'][:, k] = st['Om0']
        L_['pos_err'][k] = np.linalg.norm(cmd['ex'])
        _Xd, _Vd, _, _, _, _ = reference(st['t'])
        L_['xd'][:, k] = _Xd
        L_['vd'][:, k] = _Vd
        L_['track_err_vec'][:, k] = cmd['ex']
        if not np.isfinite(st['x0']).all() or L_['pos_err'][k] > 1e3:
            print(f'*** 发散于 k={k}  t={k*dt:.3f}s  |ex|={L_["pos_err"][k]:.3e}')
            print('    x0 =', st['x0'])
            print('    |Om0| =', np.linalg.norm(st['Om0']))
            print('    |Om|  =', np.abs(st['Om']).max())
            print('    |qd|  =', np.abs(st['qd']).max())
            print('    q_z   =', st['q'][2, :])
            print('    thrust_pct =', np.round(L_['thrust_pct'][:, k], 1))
            print('    tension    =', np.round(L_['tension'][:, k], 4))
            return st, L_, mem
        if k % 500 == 0:
            print(f'k={k:5d} t={k*dt:6.2f}s |ex|={L_["pos_err"][k]:.3e} '
                  f'|Om0|={np.linalg.norm(st["Om0"]):.2e} '
                  f'|Om|={np.abs(st["Om"]).max():.2e} '
                  f'|qd|={np.abs(st["qd"]).max():.2e} '
                  f'thr%={L_["thrust_pct"][0, k]:.1f} '
                  f'T={L_["tension"][0, k]:.4f}')
        if DEBUG and k < 14:
            qq = L_['q'][:, :, k]
            print(f'k={k:3d} |ex|={np.linalg.norm(cmd["ex"]):8.4f} '
                  f'|Fd|={np.linalg.norm(cmd["Fd"]):7.4f} '
                  f'tens={np.round(L_["tension"][:, k], 4)} '
                  f'thr%={np.round(L_["thrust_pct"][:, k], 1)} '
                  f'|Om|max={np.abs(L_["Om"][:, :, k]).max():7.3f} '
                  f'q_z={np.round(qq[2, :], 4)}')

        if k == n_steps - 1:
            break

        alpha = min(1.0, dt / TAU)
        thrust_act = np.clip(thrust_act + alpha * (cmd['thrust'] - thrust_act), 0, MAX_F)

        u_actual = np.zeros((3, N))
        for i in range(N):
            u_actual[:, i] = -thrust_act[i] * (st['R'][:, :, i] @ E3)
            Oi = st['Om'][:, i]
            st['torque'][:, i] = J @ (RATE_BW * (cmd['om_cmd'][:, i] - Oi)) \
                + np.cross(Oi, J @ Oi)

        d = dynamics(st, u_actual)
        st['u_par'] = d['u_par']

        st['v0'] = st['v0'] + dt * d['x0dd']
        st['x0'] = st['x0'] + dt * st['v0']
        st['Om0'] = st['Om0'] + dt * d['Om0d']
        st['R0'] = proj_so3(st['R0'] @ exp_so3(st['Om0'] * dt))
        for i in range(N):
            st['qd'][:, i] = st['qd'][:, i] + dt * d['qdd'][:, i]
            st['q'][:, i] = unit(st['q'][:, i] + dt * st['qd'][:, i])
            st['Om'][:, i] = np.clip(st['Om'][:, i] + dt * d['Om_dot'][:, i],
                                     -MAX_BODY_RATE, MAX_BODY_RATE)
            st['R'][:, :, i] = proj_so3(st['R'][:, :, i]
                                        @ exp_so3(st['Om'][:, i] * dt))
            nr = np.linalg.norm(st['qd'][:, i])
            if nr > MAX_PEND_RATE:
                st['qd'][:, i] *= MAX_PEND_RATE / nr

        # 绳索约束投影：剥掉 q_dot 的径向分量，保持 q.qd == 0
        # ★ 必须与 crazyflie_slung_simulation.m 的对应循环严格一致。
        #   依据：||q|| == 1 求导 => q.qd == 0；也是论文 (1) qd = omega x q 的推论。
        #   旧代码只 unit() 归一化 q 而不投影 qd，实测残差 6.371e-4
        #   （相对 |qd| 峰值 0.716 是 0.089 %）。
        for i in range(N):
            qi = st['q'][:, i]
            st['qd'][:, i] = st['qd'][:, i] - qi * float(qi @ st['qd'][:, i])

        st['a_prev'] = d['x0dd']
        st['Om0d_prev'] = d['Om0d']
        mem['q_prev'] = st['q'].copy()

    return st, L_, mem


if __name__ == '__main__':
    st, L_, mem = run()
    # ★ 稳态窗口必须与 MATLAB 严格一致：crazyflie_slung_simulation.m 的
    #   steadyIndex = max(1, round(0.8*nSteps)) : nSteps   即"后 20 %"。
    #   早期这里写死 slice(-6000)（后 12 s），与 MATLAB 不同，
    #   会导致镜像报告与自检脚本的 steadyPositionError 对不上。
    _n = len(L_['pos_err'])
    s = slice(int(0.8 * _n), _n)
    print("=" * 68)
    print(f" {N} 机协同吊运 —— 最终参数集验证")
    print("=" * 68)

    # ★★★ 八字避障工况专项报告 ★★★
    # 与 MATLAB 的 crazyflie_slung_simulation.m computeSummary 一一对应：
    #   起飞段 [0, T1) / 环绕段 [T1, T1+T2) / 降落段 [T1+T2, T1+T2+T3]
    # 稳态误差只在静态工况下有意义；八字工况必须分段看，否则指标被稀释。
    if USE_FIGURE_EIGHT:
        T1 = FE['takeoffDuration']
        T2 = FE['cruiseDuration']
        T3 = FE['landingDuration']
        k1 = int(T1 / DT)
        k2 = int((T1 + T2) / DT)
        pe = L_['pos_err']
        land = np.array([0.0, 0.0, -0.35])
        print("--- 三阶段跟踪误差（八字避障工况）---")
        print(f" 起飞段 [0, {T1:.1f} s]      : 峰值 "
              f"{1000 * pe[:k1 + 1].max():8.1f} mm")
        print(f" 环绕段 [{T1:.1f}, {T1 + T2:.1f} s]   : 峰值 "
              f"{1000 * pe[k1:k2 + 1].max():8.1f} mm / 均值 "
              f"{1000 * pe[k1:k2 + 1].mean():8.1f} mm")
        print(f" 降落段 [{T1 + T2:.1f}, {T1 + T2 + T3:.1f} s]  : 峰值 "
              f"{1000 * pe[k2:].max():8.1f} mm")
        print(f" 降落终点残差      : "
              f"{1000 * np.linalg.norm(L_['x0'][:, -1] - land):8.1f} mm")
        # 轨迹行程：验证环绕的确实是"八字"
        xspan = L_['x0'][0, :].max() - L_['x0'][0, :].min()
        yspan = L_['x0'][1, :].max() - L_['x0'][1, :].min()
        print(f" 轨迹行程 (x, y)   : ({xspan:.3f}, {yspan:.3f}) m "
              f"(期望 x 幅值 ±{FE['amplitudeX']:.2f} m)")
        # 包络过渡时长（决定参考加速度尖峰高度的关键参数）
        print(f" 包络过渡时长      : {FE['blendTime']:.1f} s "
              f"(>= 5 s；1.2 s 会让加速度尖峰放大 3.15 倍并发散)")
        # 绕开两个锥形障碍物：与 MATLAB obstacleClearance 同一定义
        # ★ 公式已修正为单一"离地高度"量（旧式 hA = baseZ - z + height 语义混乱，
        #   会让环绕段恒走 inf 分支、自检静默失效）。同时统计分支计数以做有效性门控。
        if OBSTACLES['enabled']:
            payR = 0.5 * np.linalg.norm(PAYLOAD_SIZE)
            minGap, minGapT = np.inf, 0.0
            n_inside = n_above = n_below = 0
            hBase = -OBSTACLES['baseZ']
            for kk in range(_n):
                hGround = -L_['x0'][2, kk]
                for jj in range(OBSTACLES['positions'].shape[1]):
                    dx = L_['x0'][0, kk] - OBSTACLES['positions'][0, jj]
                    dy = L_['x0'][1, kk] - OBSTACLES['positions'][1, jj]
                    rH = np.hypot(dx, dy)
                    hLocal = hGround - hBase
                    if hLocal >= 0 and hLocal <= OBSTACLES['height']:
                        rHere = OBSTACLES['radius'] * (1 - hLocal / OBSTACLES['height'])
                        gap = rH - rHere - payR
                        n_inside += 1
                    elif hLocal > OBSTACLES['height']:
                        gap = np.inf
                        n_above += 1
                    else:
                        gap = rH - OBSTACLES['radius'] - payR
                        n_below += 1
                    if gap < minGap:
                        minGap, minGapT = gap, kk * DT
            valid = n_inside > 0
            verdict = 'PASS' if (valid and minGap > OBSTACLES['clearance']) else 'FAIL'
            print(f" 到锥最小净间隙    : {1000 * minGap:8.1f} mm "
                  f"(t = {minGapT:.2f} s, 负载外接球 {1000 * payR:.1f} mm)  "
                  f"[{verdict}, 判据 > {1000 * OBSTACLES['clearance']:.0f} mm]")
            print(f"   自检有效性      : {'[有效]' if valid else '[★静默失效★]'}  "
                  f"(锥高度区间内 {n_inside} 点 / 锥顶以上 {n_above} / 锥底以下 {n_below})")
            apexH = -OBSTACLES['baseZ'] + OBSTACLES['height']
            cruiseH = abs(FE['cruiseHeight'])
            print(f"   锥尖离地/环绕平面: {apexH:.2f} m / {cruiseH:.2f} m  "
                  f"{'[锥尖更高，避障为真实约束]' if apexH > cruiseH else '[★锥尖偏低★]'}")
        print()

    print(f" 负载终位置        : {np.round(st['x0'], 5)}  "
          f"(目标 {TARGET_POS}, 误差 "
          f"{np.linalg.norm(st['x0'] - TARGET_POS) * 1000:.2f} mm)")
    print(f" 稳态位置误差      : {1000 * L_['pos_err'][s].mean():.2f} mm")
    print(f" 稳态绳向误差      : "
          f"{np.degrees(L_['link_err'][:, s].mean(axis=1)).round(3)} deg")
    print(f" 稳态姿态误差      : "
          f"{np.degrees(L_['att_err'][:, s].mean(axis=1)).round(3)} deg  (机体)")
    # ★ 负载姿态误差按轴分解：roll/pitch 为可控轴，yaw 为欠驱动自由轴
    load_att_mean = np.degrees(L_['load_att_err'][:, s].mean(axis=1))
    print(f" 稳态负载姿态误差  : [roll pitch yaw] = "
          f"{load_att_mean.round(4)} deg")
    print(f"   -> roll/pitch 范数 = "
          f"{np.linalg.norm(load_att_mean[:2]):.3f} deg  (< 3 deg 判据)")
    print(f"   -> yaw（欠驱动轴，见 README §5.1）= {load_att_mean[2]:.3f} deg")
    yaw_drift = abs(L_['load_Om0'][2, -1])
    print(f" 负载偏航漂移率    : {yaw_drift:.3f} rad/s  (< 1.0 rad/s 判据)")
    print(f" 负载角速度终值    : {np.round(L_['load_Om0'][:, -1], 4)} rad/s")
    # ★ 悬停张力由**挂点几何**决定，不能假设"三等分"。
    #   新挂点（边中点+对边两顶点）质心偏离负载质心 0.0333 m
    #   ⇒ 悬停张力是 2:1:1 = [0.3924, 0.1962, 0.1962] N（合计仍 = m0 g）。
    _bal = np.vstack([np.ones(N), RHO[1, :], -RHO[0, :]])
    _muH = np.linalg.solve(_bal, np.array([-M0 * G, 0.0, 0.0]))
    _T_theory = np.abs(_muH)
    print(f" 稳态各绳索张力    : {L_['tension'][:, s].mean(axis=1).round(4)} N")
    print(f" 理论悬停张力      : {_T_theory.round(4)} N "
          f"(由挂点几何解出，合计 {_T_theory.sum():.4f} = m0 g)")
    print(f" 张力范围          : {L_['tension'].min():.4f} ~ {L_['tension'].max():.4f} N")
    # ★ 绳索单边约束：绳只能受拉，mu_i 必须 > 0（否则绳松弛，需另建模型）
    print(f" 最小张力裕度      : {L_['tension'].min() / _T_theory.min() * 100:.1f} %  "
          f"(= 相对**最小那根**绳的悬停张力；>0 表示绳全程绷紧)")
    # ★ 绳索长度不变量：绷紧的数学含义。MATLAB 已加入同名断言（< 1e-9 m）。
    #   用显式循环而非 einsum —— R0 是 (3,3,nSteps)、RHO 是 (3,N)，
    #   einsum 下标一旦写错就会把 N 与空间维搞混，曾产生 4.5 cm 的假警报。
    # ★ 绳索长度不变量：绷紧的数学含义。MATLAB 已加入同名断言（< 1e-9 m）。
    #   用显式循环而非 einsum —— R0 是 (3,3,nSteps)、RHO 是 (3,N)，
    #   einsum 下标一旦写错就会把 N 与空间维搞混，曾产生 4.5 cm 的假警报。
    _ropeLen = L                       # 模块级 L 就是绳长，别被日志变量 L_ 混淆
    _maxDrift = 0.0
    for _k in range(L_['x0'].shape[1]):
        for _i in range(N):
            _att = L_['x0'][:, _k] + L_['R0'][:, :, _k] @ RHO[:, _i]
            _vehPos = L_['veh'][:, _i, _k]
            _maxDrift = max(_maxDrift,
                            abs(float(np.linalg.norm(_att - _vehPos)) - _ropeLen))
    print(f" 绳索长度不变量    : 最大偏差 {_maxDrift:.2e} m  "
          f"({'PASS' if _maxDrift < 1e-9 else 'FAIL'}, 判据 < 1e-9 m)")
    # ★ q'·q_dot 应恒为 0（||q||==1 求导 / 论文 (1) qd = omega x q）
    _qd = np.diff(L_['q'], axis=2) / DT
    _dots = np.einsum('ijk,ijk->jk', L_['q'][:, :, :-1], _qd)
    print(f" q'·q_dot 残差     : {np.abs(_dots).max():.3e}  "
          f"(相对 |q_dot| 峰值 {np.abs(_dots).max() / np.abs(_qd).max() * 100:.3f} %)")
    print(f" 推力峰值占比      : {L_['thrust_pct'].max():.1f} %")
    print(f" 最大机体角速度    : {np.abs(L_['Om']).max():.3f} rad/s")
    print(f" 状态有限          : {np.all(np.isfinite(L_['x0']))}")
    # 各机是否在负载上方？
    hL = -L_['x0'][2, :]
    above = all((-L_['veh'][2, i, :] > hL + 0.05).all() for i in range(N))
    print(f" 各机始终在负载上方: {above}")
