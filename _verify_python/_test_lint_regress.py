# 回归测试：确认静态检查真的能抓到当初踩过的那些"MATLAB 语义层"bug
#
# 覆盖两项：
#   ① use_before_assign —— 先用后定义（当初的 armLength）
#   ② nargout_dispatch_risk —— nargout 分派漏掉匿名句柄的 -1
#      （当初的 "此类型的变量不支持使用点进行索引"）
#
# ★ 为什么必须固化成回归测试：
#   这两类缺陷 **Python 数值镜像永远抓不到**（Python 没有"文件私有函数"、
#   没有 nargout 语义）。只有静态检查能兜住，所以检查本身也要有测试。
import os
import _lint_matlab as L

# ===========================================================================
# ① 先用后定义
# ===========================================================================

BUGGY = """function [verts, faces] = boxVerticesFaces(sizeVector)
% 局部函数私有性说明：以下 helper 各自复制一份
function out = demo(cfg)
    a = 1;
    span = 2.0;
    fprintf('span %.2f | arm %.1f %%\\n', span, 2 * armLength / span * 100);
    armLength = cfg.vehicle.armLength;
    out = a + armLength;
end
"""

OK = """function out = demo(cfg)
    armLength = cfg.vehicle.armLength;
    span = 2.0;
    fprintf('span %.2f | arm %.1f %%\\n', span, 2 * armLength / span * 100);
    out = armLength;
end
"""

print("== ① 先用后定义 ==")
print()
print("== 有 bug 的版本（应报出 armLength）==")
found = L.use_before_assign(BUGGY, set())
names = sorted({n for n, _ in found})
print("  报出:", names)
assert "armLength" in names, "回归失败：没抓到 armLength！"
print("  -> 通过")

print()
print("== 修好的版本（应干净）==")
found = L.use_before_assign(OK, set())
names = sorted({n for n, _ in found})
print("  报出:", names if names else "无")
assert "armLength" not in names, "回归失败：修好后仍误报 armLength！"
print("  -> 通过")

# ===========================================================================
# ①b catch / persistent 必须被登记为"定义"
#     （漏登记会误报；且登记必须取**最早**行号，见 register 的注释）
# ===========================================================================

PERSISTENT_CASE = """function out = demo(x)
    warnedOnce();
    out = x;
end

function warnedOnce()
    persistent fired
    if isempty(fired)
        fired = false;
    end
    if fired
        return
    end
    fired = true;
end
"""

print()
print("== ①b catch / persistent 应被登记为定义（不应误报）==")
found = L.use_before_assign(PERSISTENT_CASE, set())
names = sorted({n for n, _ in found})
print("  报出:", names if names else "无")
assert "fired" not in names, (
    "回归失败：persistent fired 被误报为先用后定义 —— "
    "说明登记没有取最早行号（正文赋值抢在了声明前面）")
print("  -> 通过")

CATCH_CASE = """function out = demo()
    try
        error('boom');
    catch err
        out = err.message;
    end
end
"""

found = L.use_before_assign(CATCH_CASE, set())
names = sorted({n for n, _ in found})
print("  catch err 报出:", names if names else "无")
assert "err" not in names, "回归失败：catch err 被误报为先用后定义"
print("  -> 通过")

# ===========================================================================
# ② nargout 分派漏掉匿名句柄的 -1
# ===========================================================================

# 当初的真实写法：只有 >= 6 / == 5 / == 4，-1 掉进 else 被当 struct。
NARGOUT_BUGGY = """function [a, b] = demo(referenceFcn, t)
numberOfOutputs = nargout(referenceFcn);
if numberOfOutputs >= 6
    [a, b, c, d, e, f] = referenceFcn(t);
elseif numberOfOutputs == 5
    [a, b, c, d, e] = referenceFcn(t);
elseif numberOfOutputs == 4
    [a, b, c, d] = referenceFcn(t);
else
    data = referenceFcn(t);
    a = data.position;
end
end
"""

# 修好的写法：显式处理 nargout < 0。
NARGOUT_FIXED = """function [a, b] = demo(referenceFcn, t)
numberOfOutputs = nargout(referenceFcn);
if numberOfOutputs < 0
    [a, b, c, d, e, f] = callReferenceAnonymous(referenceFcn, t);
    return
end
if numberOfOutputs >= 6
    [a, b, c, d, e, f] = referenceFcn(t);
    return
end
[a, b, c, d, e, f] = callReferenceSingleOutput(referenceFcn, t, []);
end
"""

print()
print("== ② nargout 分派漏掉匿名句柄的 -1 ==")
print()
print("== 有 bug 的版本（应报出分派未处理 -1）==")
risks = L.nargout_dispatch_risk(NARGOUT_BUGGY)
print("  报出:", [(fn, arg, var) for fn, arg, var, _ in risks])
assert risks, "回归失败：没抓到 nargout 分派漏掉 -1！"
assert risks[0][1] == "referenceFcn", "回归失败：抓到的变量不对"
print("  -> 通过")

print()
print("== 修好的版本（应干净）==")
risks = L.nargout_dispatch_risk(NARGOUT_FIXED)
print("  报出:", risks if risks else "无")
assert not risks, "回归失败：修好后仍误报 nargout 分派"
print("  -> 通过")

# ===========================================================================
# ②b 对照：**真实交付文件**必须也是干净的
# ===========================================================================

print()
print("== ②b 真实文件 crazyflie_slung_simulation.m 应当干净 ==")
REAL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "crazyflie_slung_simulation.m")
with open(REAL, encoding="utf-8") as fh:
    real_raw = fh.read()
real_risks = L.nargout_dispatch_risk(L.strip_comments(real_raw))
print("  nargout 分派风险:", real_risks if real_risks else "无")
assert not real_risks, "交付文件里的 nargout 分派未处理 -1！"

real_order = L.use_before_assign(real_raw, set())
real_names = sorted({n for n, _ in real_order})
print("  先用后定义:", real_names if real_names else "无")
assert not real_names, "交付文件里存在先用后定义！"
print("  -> 通过")

# ===========================================================================
# ③ 续行块语法 —— `...` 续行块里插入不以 `...` 结尾的注释行
# ===========================================================================
# ★ 由来：曾把多行注释插进 `struct(...)` 的续行块中间，MATLAB 直接报
#       "文件: crazyflie_slung_parameters.m 行: 549 列: 48  无效表达式"
#   —— 硬语法错误，却逃过了前 6 类检查（它们都不分析续行语法）。
print()
print("== ③ 续行块内的注释行必须以 ... 结尾 ==")

CONT_BAD = """cfg.vis = struct(...
    'plot', true, ...
    'animate', true, ...
    % ★ 这里插了一段注释，但没有以 ... 结尾
    % 于是语句被提前截断
    'payloadSizeScale', 1.0, ...
    'vehicleScale', 1.6);
"""
CONT_OK_OUTSIDE = """% 注释放在 struct 外面完全没问题
% 第二行也不受影响
cfg.vis = struct(...
    'plot', true, ...
    'animate', true, ...
    'payloadSizeScale', 1.0, ...
    'vehicleScale', 1.6);
"""
CONT_OK_INLINE = """cfg.vis = struct(...
    'plot', true, ...
    'animate', true);          % 行尾注释是安全的
"""


def _cont_hits(src):
    hits = []
    prev = False
    for n, l in enumerate(src.split("\n"), 1):
        s = l.strip()
        cont = L._strip_strings(l).rstrip().endswith("...")
        if prev and s.startswith("%") and not cont:
            hits.append(n)
        prev = cont
    return hits


_h_bad = _cont_hits(CONT_BAD)
print("  坏样例命中行:", _h_bad)
assert _h_bad, "续行块检查没能抓到被注释截断的语句！"

_h_o1 = _cont_hits(CONT_OK_OUTSIDE)
_h_o2 = _cont_hits(CONT_OK_INLINE)
print("  好样例（注释在 struct 外）命中行:", _h_o1)
print("  好样例（行尾注释）命中行:", _h_o2)
assert not _h_o1, "续行块检查对『注释在 struct 外』误报了！"
assert not _h_o2, "续行块检查对『行尾注释』误报了！"

# 字符串字面量里的 `...` 不能被当成续行符
_h_str = _cont_hits("a = 'x...y';\n% 跟在后面的注释\n")
print("  字符串里的 ... 不被误判:", _h_str)
assert not _h_str, "_strip_strings 没能屏蔽字符串里的 ..."

# 交付文件必须干净
for _f in ("crazyflie_slung_parameters.m", "crazyflie_slung_dynamics.m",
           "crazyflie_slung_simulation.m", "crazyflie_slung_controller.m",
           "crazyflie_slung_reference.m", "crazyflie_slung_visualization.m",
           "crazyflie_slung_demo.m", "crazyflie_slung_diagnose.m"):
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _f),
              encoding="utf-8") as fh:
        assert not _cont_hits(fh.read()), "交付文件 %s 存在续行块问题！" % _f
print("  8 个交付文件均为 0")
print("  -> 通过")

print()
print("全部回归通过：")
print("  · 先用后定义检查能抓到 armLength，且不误报已修好的代码")
print("  · catch / persistent 被正确登记为定义（且取最早行号）")
print("  · nargout 分派检查能抓到匿名句柄的 -1，且交付文件干净")
print("  · 续行块检查能抓到『注释截断语句』，且 8 个交付文件干净")
