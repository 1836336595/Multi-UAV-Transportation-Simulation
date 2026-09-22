"""锥形障碍物间隙自检的回归测试。

这个脚本的存在理由：间隙自检曾经因为**高度量的定义混乱**而静默失效 ——
所有采样点都走进了 `gap = inf` 分支，于是自检永远 PASS，但那是
"根本没检查"而不是"避障成功"。这类 bug 不报错、结果还很好看，
是最危险的一种。

本测试做两件事：
  A) 复刻"旧的错误公式"并断言它在当前锥参数下**确实恒走 inf 分支**
     —— 即证明这个坑是真实存在的，不是想象中的；
  B) 用"新的正确公式"算出真实净间隙，断言它落在预期区间内且 > 判据。

只要 B 通过而 A 复现出失效，就说明修复有效且该 bug 值得被永久看住。
"""

import numpy as np

# ---------------------------------------------------------------- 参数（与 MATLAB 一致）
A_X, A_Y = 1.50, 1.62
W_X, W_Y = 0.2 * np.pi, 0.1 * np.pi
CYCLES = 2.0
# ★ T2 18 -> 45、BLEND 6.0 -> 5.0：为让包络满幅窗口覆盖住整个叶，
#   使轨迹真正"环绕"锥（见 README §11 / _diag_orbit_period.py）
T1, T2, T3 = 3.0, 45.0, 4.0
BLEND = 5.0
CRUISE_Z = -0.55
START_X, START_Y, START_Z = 0.10, -0.06, 0.45
LAND_Z = -0.35

# 新锥参数 —— 两个锥在八字**两个叶的叶心**
OBS_POS = np.array([[-0.0187, -0.0187], [1.2252, 0.3948]])
BASE_Z = 0.00
HEIGHT = 1.05
RADIUS = 0.15
CLEARANCE_CRIT = 0.05

# 旧（错误）锥参数 —— 用来演示"旧公式 + 旧几何"如何静默失效
OLD_OBS_POS = np.array([[0.75, -0.75], [0.60, 0.60]])
OLD_BASE_Z = -0.35
OLD_HEIGHT = 0.30
OLD_RADIUS = 0.22

PAYLOAD_SIZE = np.array([0.20, 0.20, 0.020])
PAYLOAD_R = 0.5 * np.linalg.norm(PAYLOAD_SIZE)

DT = 0.002
DUR = T1 + T2 + T3


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return 6 * u**5 - 15 * u**4 + 10 * u**3


def reference_track():
    """重建参考轨迹，必须与 crazyflie_slung_reference.m / 镜像 reference() 逐式一致。

    ★★ 这里曾经有两处错误，导致重建轨迹与真实代码不符、净间隙数字偏小 ★★
       1. 包络写成了 τ 的函数 smoothstep(τ/B)*smoothstep((t_end-τ)/B)，
          而真实代码是**真实时间**的函数：blend_env(t - T1, T2, B)。
       2. 环绕段多加了 startPosition 水平偏移 (+0.10, -0.06)，
          而真实代码环绕段是 p = [p̂x·env, p̂y·env, zC]，**以原点为中心**。
    修正后本函数与镜像在 dmin 上吻合，测试才有意义。
    """
    t_end = CYCLES * 2 * np.pi / W_Y
    s = t_end / T2
    t = np.arange(0.0, DUR + DT, DT)
    px = np.zeros_like(t)
    py = np.zeros_like(t)
    pz = np.zeros_like(t)
    for k, tk in enumerate(t):
        if tk < T1:
            sg = smoothstep(tk / T1)
            px[k] = START_X + (0.0 - START_X) * sg
            py[k] = START_Y + (0.0 - START_Y) * sg
            pz[k] = START_Z + (CRUISE_Z - START_Z) * sg
        elif tk < T1 + T2:
            tau = (tk - T1) * s
            pz[k] = CRUISE_Z
            # ★ 包络是**真实时间**的函数（与 _blend_env 一致），不是 τ 的函数
            tl = tk - T1
            env = smoothstep(tl / BLEND) * smoothstep((T2 - tl) / BLEND)
            # ★ 环绕段以原点为中心，**不叠加** startPosition 水平偏移
            px[k] = env * A_X * np.sin(W_X * tau)
            py[k] = env * (A_Y / 2) * (1 - np.cos(W_Y * tau))
        else:
            sg = smoothstep((tk - T1 - T2) / T3)
            # 降落：水平回到 landPosition（本仿真是原点）
            px[k] = 0.0
            py[k] = 0.0
            pz[k] = CRUISE_Z + (LAND_Z - CRUISE_Z) * sg
    return t, px, py, pz


def old_formula_branch(px, py, pz, obs_pos=None, base_z=None, height=None,
                       radius=None):
    """复刻**旧的错误** hAboveBase 公式，返回各分支计数与 min 间隙。

    公式原文： hAboveBase = (baseZ - z) + height，再与 height 比较。
    注意它的量纲实际上是"离锥顶的距离"——baseZ 为负、height 为正时，
    z = baseZ（锥底）处给出 hAboveBase = height，z 更负（更高）时反而更大。
    所以当负载长期位于锥底以下/以下某处时它会恒 > height，全部走 inf。
    """
    obs_pos = OBS_POS if obs_pos is None else obs_pos
    base_z = BASE_Z if base_z is None else base_z
    height = HEIGHT if height is None else height
    radius = RADIUS if radius is None else radius
    n_in = n_above = n_below = 0
    mn = np.inf
    for k in range(len(pz)):
        for j in range(obs_pos.shape[1]):
            r_h = np.hypot(px[k] - obs_pos[0, j], py[k] - obs_pos[1, j])
            h_above = (base_z - pz[k]) + height          # <<< 旧公式
            if 0 <= h_above <= height:
                r_here = radius * (1 - h_above / height)
                gap = r_h - r_here - PAYLOAD_R
                n_in += 1
            elif h_above > height:
                gap = np.inf
                n_above += 1
            else:
                gap = r_h - radius - PAYLOAD_R
                n_below += 1
            mn = min(mn, gap)
    return n_in, n_above, n_below, mn


def new_formula_branch(px, py, pz, obs_pos=None, base_z=None, height=None,
                       radius=None):
    """**新的正确**公式：单一"离地高度"量。"""
    obs_pos = OBS_POS if obs_pos is None else obs_pos
    base_z = BASE_Z if base_z is None else base_z
    height = HEIGHT if height is None else height
    radius = RADIUS if radius is None else radius
    h_base = -base_z
    h_ground = -pz
    n_in = n_above = n_below = 0
    mn = np.inf
    mn_time = 0.0
    for k in range(len(pz)):
        for j in range(obs_pos.shape[1]):
            r_h = np.hypot(px[k] - obs_pos[0, j], py[k] - obs_pos[1, j])
            h_local = h_ground[k] - h_base
            if 0 <= h_local <= height:
                r_here = radius * (1 - h_local / height)
                gap = r_h - r_here - PAYLOAD_R
                n_in += 1
                if gap < mn:
                    mn = gap
                    mn_time = k * DT
            elif h_local > height:
                gap = np.inf
                n_above += 1
            else:
                gap = r_h - radius - PAYLOAD_R
                n_below += 1
            if gap < mn:
                mn = gap
    return n_in, n_above, n_below, mn, mn_time


def count_cruise_inside(px, py, pz, t, obs_pos, base_z, height, new=True):
    """只统计**环绕段** [T1, T1+T2) 内、落在锥高度区间内的采样点数。

    这是门控真正关心的量：起飞/降落段负载贴近地面，"采样到"是必然的；
    只有环绕段被检查到，才说明避障约束真的生效了。
    """
    h_base = -base_z
    n = 0
    for k in range(len(pz)):
        if not (T1 <= t[k] < T1 + T2):
            continue
        for j in range(obs_pos.shape[1]):
            if new:
                h_local = (-pz[k]) - h_base
            else:
                h_local = (base_z - pz[k]) + height     # 旧公式
            if 0 <= h_local <= height:
                n += 1
    return n


def main():
    bar = "=" * 70
    print(bar)
    print(" 锥形障碍物间隙自检 —— 回归测试")
    print(bar)
    print(" 新锥参数: 轴 (±1.50, 0.60) | 底 z=%.2f | 高 %.2f | 底半径 %.2f"
          % (BASE_Z, HEIGHT, RADIUS))
    print(" 旧锥参数: 轴 (±0.75, 0.60) | 底 z=%.2f | 高 %.2f | 底半径 %.2f"
          % (OLD_BASE_Z, OLD_HEIGHT, OLD_RADIUS))
    print(" 负载外接球半径 = %.4f m" % PAYLOAD_R)
    print()

    t, px, py, pz = reference_track()
    fail = []

    # ---------------- A) 旧的错误公式 + 旧几何：只在起飞瞬间"碰巧"检查到
    n_in_o, n_ab_o, n_be_o, mn_o = old_formula_branch(
        px, py, pz, OLD_OBS_POS, OLD_BASE_Z, OLD_HEIGHT, OLD_RADIUS)
    print("A) 旧公式 + 旧锥几何 —— 复现「只在起飞瞬间碰巧通过」:")
    print("     hAboveBase = (baseZ - z) + height = (-0.35 - z) + 0.30")
    print("     区间内 = %d | 锥顶以上(inf) = %d | 锥底以下 = %d"
          % (n_in_o, n_ab_o, n_be_o))
    print("     min 净间隙 = %s"
          % ("inf" if np.isinf(mn_o) else "%.1f mm" % (mn_o * 1000)))
    # 旧几何下的"有效采样"全部来自起飞段（负载还在地面附近），环绕段恒走 inf
    frac = n_in_o / (n_in_o + n_ab_o + n_be_o)
    print("     有效采样占比 = %.1f%%  （环绕段占全程 %.0f%%）"
          % (frac * 100, T2 / DUR * 100))
    if n_in_o > 0 and frac < 0.05:
        print("     [OK] 复现：有效采样仅 %.1f%%，全部落在起飞/降落段；"
              % (frac * 100))
        print("          环绕段（占 %.0f%%）恒走 inf 分支 ⇒ 避障从未被真正检查。"
              % (T2 / DUR * 100))
        print("          旧报告里的 560.9 mm 正是来自 t = 0.41 s 的起飞瞬间，")
        print("          与环绕段毫无关系 —— 典型的「假信心数字」。")
    else:
        fail.append("未复现旧的静默失效（有效采样占比 %.1f%%）" % (frac * 100))
        print("     [FAIL] 未复现")

    # ---------------- B) 新的正确公式 + 新锥几何
    n_in_n, n_ab_n, n_be_n, mn_n, mn_t = new_formula_branch(px, py, pz)
    print()
    print("B) 新公式 + 新锥几何 —— 真实检查:")
    print("     hLocal = 离地高度 - 锥底离地高度")
    print("     区间内 = %d | 锥顶以上(inf) = %d | 锥底以下 = %d"
          % (n_in_n, n_ab_n, n_be_n))
    print("     min 净间隙 = %.1f mm  @t = %.2f s" % (mn_n * 1000, mn_t))
    print()
    if n_in_n > 0:
        print("     [OK] 环绕段确实进入了锥的高度区间，避障被真实检查")
    else:
        fail.append("新公式下区间内采样仍为 0")
        print("     [FAIL] 区间内采样为 0")

    # 另算"只在环绕段内"的最小净间隙 —— 这才是"绕锥飞"的那个间隙。
    # ★ 全程最小间隙出现在**低空过渡段**（负载刚离地/快落地时），
    #   因为锥在底部最粗（r = 锥底半径 0.15 m），而负载当时离锥轴约 0.43 m。
    #   两者都是真实的几何约束，都要过判据，但意义不同，所以分别报。
    m_c = (t >= T1) & (t <= T1 + T2)
    n_in_c, _, _, mn_c, mn_ct = new_formula_branch(
        px[m_c], py[m_c], pz[m_c])
    print("     环绕段内 min 净间隙 = %.1f mm  (这才是绕锥飞的那个间隙)"
          % (mn_c * 1000))
    print("       全程 min 出现在 t=%.2f s，属低空过渡段（锥底部更粗）" % mn_t)

    if mn_n > CLEARANCE_CRIT:
        print("     [OK] 净间隙 %.1f mm > 判据 %.0f mm  => PASS  (%.1f 倍裕度)"
              % (mn_n * 1000, CLEARANCE_CRIT * 1000, mn_n / CLEARANCE_CRIT))
    else:
        fail.append("净间隙 %.1f mm 未达判据 %.0f mm"
                    % (mn_n * 1000, CLEARANCE_CRIT * 1000))
        print("     [FAIL] 净间隙不足")

    # 期望量级：
    #   环绕段（即"绕锥"那一段）0.3948 - 0.0714 - 0.1418 = 0.1816 m
    #   全程最小 0.4341 - 0.1500 - 0.1418 = 0.1423 m（低空过渡段，锥底最粗）
    lo, hi = 0.11, 0.22
    if lo < mn_n < hi:
        print("     [OK] 净间隙落在预期区间 (%.2f, %.2f) m 内" % (lo, hi))
    else:
        fail.append("净间隙 %.4f m 偏离预期区间 (%.2f, %.2f)" % (mn_n, lo, hi))
        print("     [FAIL] 净间隙 %.4f m 偏离预期 (%.2f, %.2f)" % (mn_n, lo, hi))
    lo_c, hi_c = 0.15, 0.22
    if lo_c < mn_c < hi_c:
        print("     [OK] 环绕段净间隙 %.1f mm 落在预期 (%.0f, %.0f) mm 内"
              % (mn_c * 1000, lo_c * 1000, hi_c * 1000))
    else:
        fail.append("环绕段净间隙 %.4f m 偏离预期 (%.2f, %.2f) m"
                    % (mn_c, lo_c, hi_c))
        print("     [FAIL] 环绕段净间隙 %.1f mm 偏离预期" % (mn_c * 1000))

    # ---------------- C) 有效性门控必须能识别"环绕段零采样"这种失效
    # 门控判据不是"全程是否采样过"，而是"环绕段是否被检查过"。
    # 旧的失效模式恰恰是：起飞段采样到了（所以全程计数不为 0），
    # 但环绕段恒走 inf 分支 —— 只看全程计数是抓不住的。
    print()
    print("C) 有效性门控检查 —— 区分「全程采样过」与「环绕段被检查过」:")
    n_in_cruise_n = count_cruise_inside(
        px, py, pz, t, OBS_POS, BASE_Z, HEIGHT, new=True)
    n_in_cruise_o = count_cruise_inside(
        px, py, pz, t, OLD_OBS_POS, OLD_BASE_Z, OLD_HEIGHT, new=False)
    print("     新公式 环绕段区间内采样 = %d" % n_in_cruise_n)
    print("     旧公式 环绕段区间内采样 = %d" % n_in_cruise_o)
    print("     新公式 全程区间内采样   = %d" % n_in_n)
    print("     旧公式 全程区间内采样   = %d  (全是起飞/降落段的)" % n_in_o)
    if n_in_cruise_n > 0 and n_in_cruise_o == 0:
        print("     [OK] 门控能区分：旧的环绕段零采样被识别，新的正常通过")
    else:
        fail.append("环绕段门控未能区分（新 %d / 旧 %d）"
                    % (n_in_cruise_n, n_in_cruise_o))
        print("     [FAIL] 门控无法区分")

    # ---------------- D) 锥尖必须高于环绕平面，避障才是真约束
    print()
    print("D) 锥尖高度 vs 环绕平面 (避障是否为真实约束):")
    apex_h = -BASE_Z + HEIGHT
    print("     锥尖离地 = -baseZ + height = %.2f + %.2f = %.2f m"
          % (-BASE_Z, HEIGHT, apex_h))
    print("     环绕平面离地 = %.2f m" % (-CRUISE_Z))
    if apex_h > -CRUISE_Z:
        print("     [OK] 锥尖高于环绕平面 %.2f m，负载必须水平绕行"
              % (apex_h - (-CRUISE_Z)))
    else:
        fail.append("锥尖未高于环绕平面，避障不是真实约束")
        print("     [FAIL] 锥尖低于环绕平面，负载会从锥顶上方飞过")

    # ---------------- E) 穿模核对：必须用**该高度处**的锥半径，而不是锥底半径
    print()
    print("E) 穿模核对（按负载所处高度处的锥半径判断，而非锥底半径）:")
    print("     负载环绕高度 = %.2f m ⇒ 该高度锥半径 r = R*(1 - h/H)" % (-CRUISE_Z))
    for j in range(OLD_OBS_POS.shape[1]):
        d = np.hypot(px - OLD_OBS_POS[0, j], py - OLD_OBS_POS[1, j]).min()
        r_here = OLD_RADIUS * (1 - (-CRUISE_Z) / OLD_HEIGHT)
        need = r_here + PAYLOAD_R
        print("     旧锥%d 轴(%+.2f,%.2f): dmin=%.4f  该高度 r=%.4f  需 > %.4f -> %s"
              % (j + 1, OLD_OBS_POS[0, j], OLD_OBS_POS[1, j], d, r_here, need,
                 "穿模" if d < need else "不穿模"))
    for j in range(OBS_POS.shape[1]):
        d = np.hypot(px - OBS_POS[0, j], py - OBS_POS[1, j]).min()
        r_here = RADIUS * (1 - (-CRUISE_Z) / HEIGHT)
        need = r_here + PAYLOAD_R
        ok = d >= need
        print("     新锥%d 轴(%+.2f,%.2f): dmin=%.4f  该高度 r=%.4f  需 > %.4f -> %s%s"
              % (j + 1, OBS_POS[0, j], OBS_POS[1, j], d, r_here, need,
                 "穿模" if not ok else "不穿模", "  [OK]" if ok else ""))
        if not ok:
            fail.append("新锥%d 几何穿模（该高度判据）" % (j + 1))
    print("     ★ 注意：若误用锥底半径 %.2f 判据（需 > %.4f）会得到「穿模」的"
          % (RADIUS, RADIUS + PAYLOAD_R))
    print("        假结论 —— 因为负载只在 0.55 m 高度飞行，永远不会靠近锥底。")
    print("        高度感知的判据才是物理正确的。")

    # ---------------- F) 环绕性核对：锥必须被轨迹**包围**，才算"围绕中心环绕"
    # ★ 这是用户新要求的核心断言，也是"绕锥飞"与"擦着锥过"的唯一区分。
    #   判定用**射线-线段求交**（不是射线-点：采样点间距约 50 mm，
    #   用"垂距 < 4 mm"当命中会恒判未包围 —— 这个测试自身的 bug 曾误导过我）。
    print()
    print("F) 环绕性核对（锥是否真的被轨迹包围）：")
    m_cruise = (t >= T1) & (t <= T1 + T2)
    cpath = np.column_stack([px[m_cruise], py[m_cruise]])

    def surrounded(cxy, pts, nray=32):
        A = pts[:-1]
        E = pts[1:] - pts[:-1]
        p = A - cxy
        for ang in np.linspace(0, 2 * np.pi, nray, endpoint=False):
            d = np.array([np.cos(ang), np.sin(ang)])
            den = d[0] * E[:, 1] - d[1] * E[:, 0]
            ok = np.abs(den) > 1e-14
            tt = np.full(len(E), -1.0)
            ss = np.full(len(E), -1.0)
            tt[ok] = (p[ok, 0] * E[ok, 1] - p[ok, 1] * E[ok, 0]) / den[ok]
            ss[ok] = (p[ok, 0] * d[1] - p[ok, 1] * d[0]) / den[ok]
            if not np.any((tt > 1e-6) & (ss >= 0) & (ss <= 1)):
                return False
        return True

    for j in range(OBS_POS.shape[1]):
        cxy = OBS_POS[:, j]
        ok = surrounded(cxy, cpath)
        d = np.hypot(px - cxy[0], py - cxy[1]).min()
        print("     锥%d 轴(%+.4f,%+.4f): 被轨迹包围 = %s | 到轨迹最近 %.4f m"
              % (j + 1, cxy[0], cxy[1], "是" if ok else "否", d))
        if not ok:
            fail.append("锥%d 未被轨迹包围 —— 只是擦过而非环绕" % (j + 1))
    if not fail or all("包围" not in f for f in fail):
        print("     [OK] 两个锥都被轨迹包围 ⇒ 负载确实在**绕锥飞**")

    print()
    print(bar)
    if fail:
        print(" 结果: FAIL")
        for f in fail:
            print("   - " + f)
    else:
        print(" 结果: 全部通过 —— 静默失效已复现，新公式给出真实净间隙且 PASS")
    print(bar)
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
