# 验证 callReferenceFunction 的分派逻辑（逻辑仿真，不是数值仿真）
#
# 背景：本机不能运行 MATLAB，而这次的缺陷恰恰只在 MATLAB 语义层暴露 ——
#       **匿名函数句柄的 nargout() 恒返回 -1**，旧分派只写 >=6 / ==5 / ==4，
#       于是 -1 掉进最后的 struct 分支，报
#       "此类型的变量不支持使用点进行索引"。
#
# 本脚本做两件事：
#   A) 结构审计：从真实 .m 文件里确认修复的结构确实存在（-1 守卫在前、
#      匿名路径先试 6、兜底走 struct 判定）。
#   B) 逻辑仿真：用一个"会按 MATLAB 规则抛 Too-many-output 的假参考函数"
#      驱动一份分派逻辑镜像，对 7 种契约 × 各种 nargout 取值断言输出正确。
#
# ★ 明确边界：B 验证的是**分派逻辑本身**（回退次序、默认值补齐、不出现未赋值
#   输出），它并不能替代真跑 MATLAB。它所依赖并固化的 MATLAB 语义假设是：
#   **body 为单个函数调用的匿名函数会转发全部输出**，即
#       f = @(t) g(t);
#       [a,b,c,d,e,f] = f(t)      % 等价于 g 返回 6 个
#   这正是本项目的 cfg.referenceFcn = @(t) crazyflie_slung_reference(t, cfg)。
import os
import re

import _lint_matlab as L

ZERO3 = (0.0, 0.0, 0.0)
EYE3 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

HERE = os.path.dirname(os.path.abspath(__file__))
REAL_M = os.path.join(HERE, "..", "crazyflie_slung_simulation.m")


class TooManyOutputs(Exception):
    """MATLAB 的 "Too many output arguments" —— 请求的输出多于函数提供的。"""


class FakeReferenceFcn:
    """按 MATLAB 规则响应输出个数请求的假参考函数。"""

    def __init__(self, n_max, struct_mode=False, struct_fields=None):
        self.n_max = n_max
        self.struct_mode = struct_mode
        self.struct_fields = struct_fields or set()
        self.calls = []

    def call(self, nout):
        self.calls.append(nout)
        if self.struct_mode:
            if nout > 1:
                raise TooManyOutputs()
            return {"struct": True}
        if nout > self.n_max:
            raise TooManyOutputs()
        vals = [
            (1.0, 2.0, 3.0),                       # position
            (0.1, 0.2, 0.3),                       # velocity
            (0.9, 0.8, 0.7),                       # acceleration
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),  # rotation
            (0.01, 0.02, 0.03),                    # bodyRate
            (0.05, 0.06, 0.07),                    # bodyRateDot
        ]
        return tuple(vals[:nout])


# ---------------------------------------------------------------------------
# 分派逻辑的镜像（与 crazyflie_slung_simulation.m 的 callReferenceFunction 一一对应）
# ---------------------------------------------------------------------------

class _Log:
    def __init__(self):
        self.degraded = []


def call_reference_function(fn, nargout_value, log):
    if nargout_value < 0:
        return _call_anonymous(fn, log)

    body_rate = ZERO3
    body_rate_dot = ZERO3
    if nargout_value >= 6:
        r = fn.call(6)
        return r[0], r[1], r[2], r[3], r[4], r[5]
    if nargout_value == 5:
        log.degraded.append("5")
        r = fn.call(5)
        return r[0], r[1], r[2], r[3], r[4], body_rate_dot
    if nargout_value == 4:
        log.degraded.append("4")
        r = fn.call(4)
        return r[0], r[1], r[2], r[3], body_rate, body_rate_dot
    if nargout_value == 3:
        log.degraded.append("3")
        r = fn.call(3)
        return r[0], r[1], r[2], EYE3, body_rate, body_rate_dot
    return _call_single_output(fn, log, first_error=None)


def _call_anonymous(fn, log):
    body_rate = ZERO3
    body_rate_dot = ZERO3
    first_error = None
    try:
        r = fn.call(6)
        return r[0], r[1], r[2], r[3], r[4], r[5]
    except TooManyOutputs as exc:
        first_error = exc
    try:
        log.degraded.append("5")
        r = fn.call(5)
        return r[0], r[1], r[2], r[3], r[4], body_rate_dot
    except TooManyOutputs:
        log.degraded.pop()
    try:
        log.degraded.append("4")
        r = fn.call(4)
        return r[0], r[1], r[2], r[3], body_rate, body_rate_dot
    except TooManyOutputs:
        log.degraded.pop()
    try:
        log.degraded.append("3")
        r = fn.call(3)
        return r[0], r[1], r[2], EYE3, body_rate, body_rate_dot
    except TooManyOutputs:
        log.degraded.pop()
    return _call_single_output(fn, log, first_error)


def _call_single_output(fn, log, first_error):
    body_rate = ZERO3
    body_rate_dot = ZERO3
    try:
        data = fn.call(1)
    except TooManyOutputs:
        raise first_error if first_error is not None else RuntimeError("uncallable")
    # MATLAB 里"单个输出"就是裸值，不是 1-元组；假函数为了统一返回元组，
    # 这里按 MATLAB 语义拆开。
    if isinstance(data, tuple) and len(data) == 1:
        data = data[0]
    if isinstance(data, dict) and data.get("struct"):
        position = (1.0, 2.0, 3.0)
        velocity = (0.1, 0.2, 0.3)
        acceleration = (0.9, 0.8, 0.7)
        rotation = EYE3
        if not {"bodyRate", "bodyRateDot"} <= fn.struct_fields:
            log.degraded.append("struct")
        return position, velocity, acceleration, rotation, body_rate, body_rate_dot
    # ★ 与 MATLAB 源码的 `else` 分支逐行对应：退化到位置向量时**必须告警**
    log.degraded.append("1")
    return data, ZERO3, ZERO3, EYE3, body_rate, body_rate_dot


# ---------------------------------------------------------------------------
# A) 结构审计：真实文件里修复必须存在
# ---------------------------------------------------------------------------

print("=" * 74)
print(" A) 结构审计 —— 真实 crazyflie_slung_simulation.m")
print("=" * 74)

with open(REAL_M, encoding="utf-8") as fh:
    src_raw = fh.read()
# ★ 用项目自己的 strip_code 去注释与字符串，而不是 `%[^\n]*` 粗剥 ——
#   后者会把 warning 格式串里的 `%s` 之后整行吃掉，破坏代码结构。
src = L.strip_code(src_raw)

checks = []

# A1 必须有显式的 nargout < 0 守卫，且出现在 >= 6 之前
i_neg = src.find("if numberOfOutputs < 0")
i_ge6 = src.find("if numberOfOutputs >= 6")
checks.append(("存在 `if numberOfOutputs < 0` 守卫", i_neg != -1))
checks.append(("`< 0` 守卫出现在 `>= 6` 之前", i_neg != -1 and i_ge6 != -1 and i_neg < i_ge6))

# A2 匿名路径的 helper 必须存在并被调用
checks.append(("定义了 callReferenceAnonymous",
               re.search(r"function\b[^\n]*callReferenceAnonymous\s*\(", src) is not None))
checks.append(("`< 0` 分支调用了 callReferenceAnonymous",
               re.search(r"if numberOfOutputs < 0.*?callReferenceAnonymous", src, re.S) is not None))

# A3 匿名路径里必须逐级回退尝试 6 / 5 / 4 / 3 个输出
#    ★ 两个坑：
#      (1) strip_code 会把 `...` 续行合并成一行，所以不能按"函数头换行"匹配；
#      (2) 边界不能写 `.*?\nend\n` —— 匿名函数体内部的 try/catch 也有 `end`，
#          非贪婪会在第一个 catch 的 end 处截断（实测踩过）。
#      正解：定位到函数头，再切到**下一个 function 头**为止。
mh = re.search(r"function\b[^\n]*callReferenceAnonymous\s*\([^)]*\)", src)
anon_body = ""
if mh:
    rest = src[mh.end():]
    nxt = rest.find("\nfunction ")
    anon_body = rest[:nxt] if nxt != -1 else rest
n_try_calls = anon_body.count("referenceFcn(t)")
checks.append((f"匿名路径逐级回退（实际尝试 {n_try_calls} 次，应 >= 4）", n_try_calls >= 4))
checks.append(("匿名路径第一发就请求 6 个输出",
               re.search(r"\[position, velocity, acceleration, rotation, bodyRate, bodyRateDot\]"
                         r"\s*=\s*referenceFcn\(t\)", anon_body) is not None))

# A4 必须有单输出兜底与 struct 判定
checks.append(("定义了 callReferenceSingleOutput", "callReferenceSingleOutput(referenceFcn, t," in src_raw))
checks.append(("单输出兜底里有 isstruct 判定",
               re.search(r"callReferenceSingleOutput\(referenceFcn, t, firstError\).*?isstruct\(data\)", src, re.S) is not None))

# A5 退化时必须告警（不能静默）
checks.append(("有 warnDegraded 一次性告警", "function warnDegraded(kind)" in src_raw))

# A6 真实句柄确实是匿名函数（这正是 -1 的来源）
with open(os.path.join(HERE, "..", "crazyflie_slung_parameters.m"), encoding="utf-8") as fh:
    p_raw = fh.read()
anon_handle = re.search(r"cfg\.referenceFcn\s*=\s*@\s*\(", p_raw)
checks.append(("cfg.referenceFcn 是匿名句柄（nargout 必为 -1）", anon_handle is not None))

ok_all = True
for name, ok in checks:
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}")
    ok_all = ok_all and ok

assert ok_all, "结构审计失败：修复不完整"

# ---------------------------------------------------------------------------
# B) 逻辑仿真：7 种契约 × 各 nargout 取值
# ---------------------------------------------------------------------------

print()
print("=" * 74)
print(" B) 逻辑仿真 —— 分派在每种契约下都必须给出完整、正确的输出")
print("=" * 74)

cases = [
    ("6 输出（本项目实际契约）", FakeReferenceFcn(6), -1, (0.05, 0.06, 0.07), []),
    ("5 输出", FakeReferenceFcn(5), -1, ZERO3, ["5"]),
    ("4 输出", FakeReferenceFcn(4), -1, ZERO3, ["4"]),
    ("3 输出", FakeReferenceFcn(3), -1, ZERO3, ["3"]),
    ("单输出 struct（含两个速率字段）",
     FakeReferenceFcn(1, struct_mode=True, struct_fields={"bodyRate", "bodyRateDot"}),
     -1, ZERO3, []),
    ("单输出 struct（缺两个速率字段）",
     FakeReferenceFcn(1, struct_mode=True, struct_fields=set()),
     -1, ZERO3, ["struct"]),
    ("单输出位置向量",
     FakeReferenceFcn(1), -1, ZERO3, ["1"]),
]

failures = []
for name, fn, nout, want_rate_dot, want_degraded in cases:
    log = _Log()
    pos, vel, acc, rot, rate, rate_dot = call_reference_function(fn, nout, log)
    problems = []
    for label, val in (("position", pos), ("velocity", vel), ("acceleration", acc),
                       ("rotation", rot), ("bodyRate", rate), ("bodyRateDot", rate_dot)):
        if val is None:
            problems.append(f"{label} 未赋值")
    if rate_dot != want_rate_dot:
        problems.append(f"bodyRateDot 期望 {want_rate_dot} 实得 {rate_dot}")
    if log.degraded != want_degraded:
        problems.append(f"退化告警期望 {want_degraded} 实得 {log.degraded}")
    status = "OK " if not problems else "FAIL"
    print(f"  [{status}] {name:34s} nargout={nout:2d}  "
          f"首次请求输出数={fn.calls[0]}  告警={log.degraded}")
    if problems:
        for pr in problems:
            print(f"         - {pr}")
        failures.append(name)

# 关键断言：本项目真实契约下**必须第一发命中 6 输出**，不得退化
fn6 = FakeReferenceFcn(6)
log6 = _Log()
call_reference_function(fn6, -1, log6)
print()
print(f"  ★ 本项目实际契约：nargout=-1 -> 首次请求 {fn6.calls[0]} 个输出，"
      f"共调用 {len(fn6.calls)} 次，退化告警 {log6.degraded}")
assert fn6.calls[0] == 6, "匿名句柄路径没有优先请求 6 个输出！"
assert len(fn6.calls) == 1, "匿名句柄路径下 6 输出契约不应回退试探"
assert log6.degraded == [], "6 输出契约不应产生退化告警"
print("  -> 第一发命中，无回退、无告警")

# 真实契约下 nargout 已知的分支也应干净
for nv in (6, 7):
    log = _Log()
    fn = FakeReferenceFcn(6)
    call_reference_function(fn, nv, log)
    assert log.degraded == [], f"nargout={nv} 不应告警"
    assert fn.calls == [6], f"nargout={nv} 应请求 6 输出"
print("  -> nargout=6/7 走精确分派，同样干净")

assert not failures, f"逻辑仿真失败: {failures}"

print()
print("=" * 74)
print(" 全部通过：")
print("   · 修复结构完整（-1 守卫在前、匿名逐级回退、struct 兜底、退化告警）")
print("   · 7 种契约在 nargout=-1 下都返回完整输出，默认值补齐符合预期")
print("   · 本项目真实契约第一发命中 6 输出，不回退也不告警")
print("=" * 74)
