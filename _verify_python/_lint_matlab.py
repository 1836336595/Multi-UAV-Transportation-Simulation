"""MATLAB static lint for the v3 slung-load code.

Catches the class of bug that a Python numerical mirror CANNOT catch:
  * calling a helper that is not defined in the SAME file (MATLAB local
    functions are file-private -> cross-file calls fail at runtime)
  * dotted field paths (cfg.a.b / state.a / command.a) that do not exist
  * identifiers used but never assigned

Usage:  python _lint_matlab.py
"""
import re
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.abspath(os.path.join(ROOT, ".."))

FILES = [
    "crazyflie_slung_parameters.m",
    "crazyflie_slung_reference.m",
    "crazyflie_slung_dynamics.m",
    "crazyflie_slung_controller.m",
    "crazyflie_slung_simulation.m",
    "crazyflie_slung_visualization.m",
    "crazyflie_slung_demo.m",
    "crazyflie_slung_diagnose.m",
]
# ★★ 血的教训：`crazyflie_slung_reference.m` 一度**不在 FILES 里**，只被登记成
#    "同目录兄弟文件"用于可见性检查，于是成了唯一的静态检查盲区 ——
#    而它恰好是**所有参考轨迹参数（aY / cruiseDuration / blendTime / 三阶段时长）
#    唯一真正落地的地方**。该文件里任何"先用后定义 / 块配平 / 名字拼写"问题
#    都逃过了检查。现已纳入逐文件检查。
SIBLING_M_FILES = [
    "crazyflie_slung_reference.m",
]

KEYWORDS = {
    "if", "elseif", "else", "end", "for", "while", "switch", "case",
    "otherwise", "try", "catch", "function", "return", "break", "continue",
    "global", "persistent", "parfor", "spmd", "classdef", "properties",
    "methods", "events", "enumeration", "arguments", "true", "false",
    # MATLAB 命令式语法里的"模式字"：`axis equal` / `yyaxis left` /
    # `grid on` / `hold on` / `view(...)` 的参数。它们不是变量，若不加进来
    # 会被"先用后定义"检查误报（实测踩过）。
    "equal", "tight", "square", "normal", "vis3d", "manual", "auto", "image",
    "left", "right", "on", "off", "xy", "x", "y", "z", "ij", "filled",
    "flat", "interp", "phong", "point", "first", "last", "none",
}

# base-MATLAB functions (R2021b) used or plausibly used by this project.
BUILTINS = set("""
abs acos acosh acot acoth acsc acsch angle asec asech asin asinh atan atan2
atanh ceil cos cosh cot coth csc csch exp fix floor hypot log log10 log2
mod power rem round sec sech sign sin sinh sqrt tan tanh nthroot
colon cumprod cumsum diff dot cross kron max min prod sort sum
length ndims numel size isempty isequal isequaln isa isfield isstruct
isfinite isinf isnan isreal isrow iscolumn isscalar isvector ismatrix
zeros ones eye nan inf eps pi true false rand randn randi
reshape repmat cat horzcat vertcat permute ipermute squeeze shiftdim
transpose ctranspose inv det trace rank svd eig eigvals qr lu chol pinv
norm null orth expm logm sqrtm mldivide mrdivide
linspace logspace meshgrid ndgrid interp1 interp2
deal cell struct struct2cell cell2mat mat2cell num2cell fieldnames.
char string strcmp strcmpi strcat strjoin strsplit sprintf sscanf regexp
regexprep strfind contains startsWith endsWith pad num2str str2num str2double
double single int8 int16 int32 int64 uint8 uint16 uint32 uint64 logical cast
disp display fprintf error warning assert nargin nargout narginchk nargoutchk
inputName nargin isargout validateattributes
rethrow throw lastwarn errordlg warndlg lasterr lasterror
tic toc cputime pause clock etime now datestr datenum datetime seconds
figure hold plot plot3 stem plotmatrix grid axis xlabel ylabel zlabel
title legend subplot xlim ylim zlim view drawnow refresh getframe
line patch surface text annotation colorbar colormap shading lighting
bar barh histogram scatter scatter3 quiver quiver3 contour contour3
fill area semilogy semilogx loglog errorbar box on off daspect pbaspect
campos camtarget camva camup camlight camproj cla clf close gca gcf gco
axes yyaxis linkaxes
gobjects VideoWriter writeVideo open waitbar
yline xline xticks yticks zticks xticklabels yticklabels
movegui uicontrol uitable uipanel warndlg errordlg msgbox
zeros ones cell repmat reshape
deg2rad rad2deg cart2sph sph2cart cart2pol pol2cart
all any find flip fliplr flipud rot90 triu tril diag blkdiag
conv deconv filter fft ifft interp
cumtrapz trapz integral quad gauss
poly polyfit polyval roots conv
optimset fminsearch fzero fminbnd ode45 ode15s odeset
who whos clear exist which path addpath rmpath cd pwd dir ls
clc format num2str eval feval builtin nargin
copyfile mkdir
syms assume simplify solve subs
plotRotation
istable table array2table readtable writetable
jsonencode jsondecode
strtrim strrep upper lower
normalize vecnorm
fieldnames real imag mean median std var set get getfield setfield rmfield
orderfields structfun cellfun arrayfun isletter isspace isnumeric ischar
isstring islogical iscell isobject ishandle strlength nnz issorted unique ind2sub fileread fullfile fclose fopen class tempdir
cummax cummin maxk mink bounds trapz gradient movmean movmedian
cellstr num2cell cell2struct struct2cell
rcond mat2str norm rank det trace linspace repmat arrayfun
""".split())

# Identifiers that are legitimately called as variables holding function handles
CALLABLE_VARS = {
    "referenceFcn", "cfg.referenceFcn", "userCfg", "options",
}


def strip_comments(text):
    """Remove only comments (keep string literals) -> for field-name extraction."""
    text = re.sub(r"\.\.\.\s*\n", " ", text)
    out = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "%":
            j = text.find("\n", i)
            if j == -1:
                j = n
            out.append(" " * (j - i))
            i = j
        elif ch == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
        elif ch == "'":
            k = i - 1
            while k >= 0 and text[k] in " \t":
                k -= 1
            if k >= 0 and (text[k].isalnum() or text[k] in "_)]}."):
                out.append("'")
                i += 1
            else:
                j = i + 1
                while j < n:
                    if text[j] == "'":
                        if j + 1 < n and text[j + 1] == "'":
                            j += 2
                            continue
                        break
                    if text[j] == "\n":
                        break
                    j += 1
                out.append(text[i:j + 1])
                i = j + 1
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def strip_code(text):
    """Remove comments AND string literals; keep single-quote transpose."""
    # join line continuations first
    text = re.sub(r"\.\.\.\s*\n", " ", text)
    out = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "%":
            # comment to end of line
            j = text.find("\n", i)
            if j == -1:
                j = n
            out.append(" " * (j - i))
            i = j
        elif ch == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 1
            out.append(" " * (j - i + 1))
            i = j + 1
        elif ch == "'":
            # transpose if previous non-space char is alnum/_/)/]/}/.
            k = i - 1
            while k >= 0 and text[k] in " \t":
                k -= 1
            if k >= 0 and (text[k].isalnum() or text[k] in "_)]}."):
                out.append("'")          # transpose, keep
                i += 1
            else:
                # string literal (doubled '' escapes)
                j = i + 1
                while j < n:
                    if text[j] == "'":
                        if j + 1 < n and text[j + 1] == "'":
                            j += 2
                            continue
                        break
                    if text[j] == "\n":      # unterminated string -> be safe
                        break
                    j += 1
                out.append(" " * (j - i + 1))
                i = j + 1
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def nargout_dispatch_risk(nocmt_code):
    """Find `nargout(x)` / `nargin(x)` dispatches that never handle the -1 case.

    ★ 为什么单独立一项检查：
      MATLAB 对**匿名函数句柄**执行 nargout() 恒返回 -1（输出个数"未知"）。
      而本项目 cfg.referenceFcn = @(t) crazyflie_slung_reference(t, cfg) 正是
      匿名句柄。旧版分派写成 `>= 6 / == 5 / == 4 / else`，-1 三个都不命中 →
      掉进最后的 else 被当成 struct，报"此类型的变量不支持使用点进行索引"。
      这类缺陷 Python 数值镜像**永远抓不到**（Python 没有 nargout 语义），
      只有静态检查能兜住，所以必须固化成一条永久规则。

    判定：对每个 `var = nargout(<标识符>)` / `var = nargin(<标识符>)`，
    要求同文件内出现对 `var` 的负数守卫（`var < 0` / `var <= -1` / `var == -1`）。
    `nargout(@name)` 这类**字面句柄**不受影响（命名函数的 nargout 已知），
    因为正则要求括号内是裸标识符。
    """
    risks = []
    for m in re.finditer(r"(\w+)\s*=\s*(nargout|nargin)\s*\(\s*([A-Za-z_]\w*)\s*\)",
                         nocmt_code):
        var, fn, arg = m.group(1), m.group(2), m.group(3)
        handled = re.search(
            r"\b%s\s*(?:<\s*0\b|<=?\s*-\s*1\b|==\s*-\s*1\b)" % re.escape(var),
            nocmt_code)
        if not handled:
            line = nocmt_code[:m.start()].count("\n") + 1
            risks.append((fn, arg, var, line))
    return risks


def use_before_assign(raw_code, local_funs):
    """Find variables read before any assignment in the same function scope.

    This is the check that would have caught `armLength` being used in a
    fprintf before its `armLength = ...` line: MATLAB raises
    "函数或变量 'armLength' 无法识别" only at run time, and the numeric
    Python mirror never sees it because the mirror is written independently.

    ★ 前提：必须传入**已剥离字符串与注释**的代码（strip_code 的输出）。
    否则中文字符串里的 "pitch" / "yaw"、属性名 'Location' 等都会被当成变量，
    产生大量假阳性。实测早期版本就踩了这个坑。
    代码内部的 strip_code 会合并 `...` 续行，因此行号是"合并后"的行号，
    与源文件略有偏移 —— 对本检查（只看先后次序）没有影响。
    """
    code = strip_code(raw_code)
    lines = code.split("\n")
    # locate function bodies: (start_line, name) for each `function ...` header
    headers = []
    for idx, ln in enumerate(lines):
        m = re.match(r"\s*function\b.*", ln)
        if m:
            nm = None
            mm = re.search(r"=\s*(\w+)\s*\(", ln) or re.search(r"function\s+(\w+)\s*\(", ln)
            if mm:
                nm = mm.group(1)
            headers.append((idx, nm))

    def scope_of(idx):
        owner = None
        for start, nm in headers:
            if start <= idx:
                owner = (start, nm)
            else:
                break
        return owner

    # first assignment line per (scope, name)
    # ★ 必须能识别**行内**的赋值（`roll = rpy(1); pitch = rpy(2);` 里的 pitch），
    #   所以不能只做 ^ 锚定的整行匹配。
    # ★ 切分时必须保护方括号内部：`[U, ~, V] = svd(R)` 若按逗号切，
    #   会碎成 `[U` / ` ~` / ` V] = svd(R)`，丢掉了 U 和 V（实测踩过）。
    #   做法：先按顶层分号切，再对每段取第一个顶层 `=` 的左侧。
    first_assign = {}

    def register(scope_and_tok, line):
        """登记某个 (scope, name) 的**最早**定义行。

        ★ 不能用 setdefault：登记是分多个循环进行的，正文赋值（`firedKinds = {}`）
          与声明（`persistent firedKinds`）散落在不同循环里。setdefault 会让
          "先跑到的循环"胜出而不是"行号最小的胜出"，于是声明行 386 被记成 388，
          反过来把 387 行的首次使用误报成"先用后定义"（实测踩过）。
          取 min 才是"首次定义"的正确语义。
        """
        prev = first_assign.get(scope_and_tok)
        if prev is None or line < prev:
            first_assign[scope_and_tok] = line

    def top_level_assignments(ln):
        """Return the LHS variable names of every `=` at bracket depth 0."""
        depth = 0
        seg_start = 0
        out = []
        i = 0
        while i <= len(ln):
            ch = ln[i] if i < len(ln) else ";"
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == ";" and depth == 0:
                seg = ln[seg_start:i]
                seg_start = i + 1
                out.append(seg)
            i += 1
        lhs_names = []
        for seg in out:
            # 在深度 0 找第一个 '='（排除 ==, <=, >=, ~=）
            d2 = 0
            for j, c in enumerate(seg):
                if c in "([{":
                    d2 += 1
                elif c in ")]}":
                    d2 -= 1
                elif c == "=" and d2 == 0:
                    if j + 1 < len(seg) and seg[j + 1] == "=":
                        break
                    if j > 0 and seg[j - 1] in "<>~=":
                        break
                    lhs = seg[:j]
                    for tok in re.findall(r"[A-Za-z_]\w*", lhs):
                        lhs_names.append(tok)
                    break
        return lhs_names

    for idx, ln in enumerate(lines):
        scope = scope_of(idx)
        body = re.sub(r"^\s*function\b", "", ln)
        for tok in top_level_assignments(body):
            register((scope, tok), idx)
    # loop vars and struct-field bases are assignments too
    for idx, ln in enumerate(lines):
        scope = scope_of(idx)
        for tok in re.findall(r"\bfor\s+(\w+)\s*=", ln):
            register((scope, tok), idx)
        for tok in re.findall(r"\b(\w+)\.\w+\s*=", ln):
            register((scope, tok), idx)
        # ★ 匿名函数参数 `@(t) ...` / `@(x, y) ...` 是**绑定变量**，
        #   在函数体内部属于已定义。不登记它们会把 `@(t) cfg.f(t)` 里的 t
        #   报成"先用后定义"（实测在 parameters.m:361 与 demo.m:58 误报）。
        #   登记位置用本行的行号即可，因为参数的作用域就是同一个表达式。
        for grp in re.findall(r"@\s*\(([^)]*)\)", ln):
            for tok in grp.split(","):
                tok = tok.strip()
                if re.fullmatch(r"\w+", tok):
                    register((scope, tok), idx)
        # ★ `catch err` 里的 err 是**绑定变量**（catch 子句即其定义处）。
        #   不登记会把 catch 块内首次使用 err 报成"先用后定义"（实测踩过）。
        for tok in re.findall(r"\bcatch\s+(\w+)", ln):
            register((scope, tok), idx)
        # ★ `persistent x y` / `global x y` 同样是定义（可空格或逗号分隔多个）。
        #   实测漏登记会把 `persistent firedKinds` 报成先用后定义。
        for kw in ("persistent", "global"):
            for m2 in re.finditer(r"\b%s\s+([^;=]*)" % kw, ln):
                for tok in re.split(r"[,\s]+", m2.group(1).strip()):
                    if re.fullmatch(r"\w+", tok):
                        register((scope, tok), idx)

    # function arguments are assigned at the header
    for start, _ in headers:
        m = re.search(r"\(([^)]*)\)", lines[start])
        if m:
            scope = scope_of(start)
            for tok in m.group(1).split(","):
                tok = tok.strip()
                if re.fullmatch(r"\w+", tok):
                    register((scope, tok), start)
        # return vars of `function [a,b] = f`
        mr = re.match(r"\s*function\s*\[([^\]]*)\]", lines[start])
        if mr:
            scope = scope_of(start)
            for tok in mr.group(1).split(","):
                tok = tok.strip()
                if re.fullmatch(r"\w+", tok):
                    register((scope, tok), start)
        mr1 = re.match(r"\s*function\s+(\w+)\s*=", lines[start])
        if mr1:
            scope = scope_of(start)
            register((scope, mr1.group(1)), start)

    findings = []
    seen = set()
    # 词法：完整的标识符（不能截断！早期版本用 [A-Za-z_]\w* 配合错误的
    # 否定后行断言，把 `crazyflie_slung_parameters` 截成 `crazyflie_slung_parameter`，
    # 制造了一堆假阳性）。这里用严格的边界 + 排除紧邻 ' 或 . 的情况。
    ident_re = re.compile(r"(?<![\w.'])([A-Za-z_]\w*)(?![\w.])")
    # 字符串字面量（MATLAB 单引号）。注意 MATLAB 里 ' 也是转置运算符，
    # 这里只做"行内成对字符串"的粗剥离 —— 用于屏蔽中文说明与 'Location' 之类的
    # 属性名，避免它们被当成变量。
    str_re = re.compile(r"'[^']*'")
    for idx, ln in enumerate(lines):
        scope = scope_of(idx)
        # 函数头行：参数已在 first_assign 里登记，且形如 function f(a,b,cfg)
        # 的行本身不构成"使用"，跳过。
        if re.match(r"\s*function\b", ln):
            continue
        stripped = str_re.sub(" ", ln)          # 抹掉字符串，保留结构
        # 抹掉 '...' 续行符残留
        stripped = stripped.replace("...", " ")
        for m in ident_re.finditer(stripped):
            tok = m.group(1)
            if tok in KEYWORDS or tok in BUILTINS or tok in local_funs:
                continue
            if tok in CALLABLE_VARS:
                continue
            # ★ 紧跟 '.' 的是结构体字段名，不是变量（`cfg.position` 里的
            #   position）。ident_re 的后向断言已排除前一字符是 '.' 的情况，
            #   这里再防一手"字段名后面紧跟 ."的链式访问。
            if m.end() < len(stripped) and stripped[m.end()] == ".":
                continue
            # 紧跟 '(' 的是函数调用，不参与变量检查
            if m.end() < len(stripped) and stripped[m.end()] == "(":
                continue
            key = (scope, tok)
            fa = first_assign.get(key)
            if fa is None or fa > idx:
                if key in seen:
                    continue
                seen.add(key)
                findings.append((tok, idx + 1))
    return findings


def parse_file(path):
    raw = open(path, encoding="utf-8").read()
    code = strip_code(raw)
    nocmt = strip_comments(raw)

    local_funs = set(re.findall(r"^\s*function\s+(\w+)", code, re.M))
    # also `function [a,b] = name(...)` and `function a = name(...)`
    for m in re.finditer(r"^\s*function\s+(?:\[[^\]]*\]|\w+)\s*=\s*(\w+)", code, re.M):
        local_funs.add(m.group(1))
    for m in re.finditer(r"^\s*function\s+(\w+)\s*\(", code, re.M):
        local_funs.add(m.group(1))

    # assigned variables: LHS of `=` (first token) and function signature args
    assigned = set()
    for m in re.finditer(r"^\s*(?:function\s+)?(?:\[([^\]]*)\]|(\w+))\s*=", code, re.M):
        for grp in m.groups():
            if grp:
                for tok in grp.split(","):
                    tok = tok.strip()
                    if re.fullmatch(r"\w+", tok):
                        assigned.add(tok)
    # function input args
    for m in re.finditer(r"^\s*function\s+(?:\[[^\]]*\]|\w+)\s*=\s*\w+\s*\(([^)]*)\)", code, re.M):
        for tok in m.group(1).split(","):
            tok = tok.strip()
            if re.fullmatch(r"\w+", tok):
                assigned.add(tok)
    # loop variables
    assigned |= set(re.findall(r"\bfor\s+(\w+)\s*=", code))
    # struct outputs assigned per field, e.g. cfg.x = ...
    assigned |= set(re.findall(r"\b(\w+)\.\w+\s*=", code))

    calls = defaultdict(list)
    for m in re.finditer(r"(?<![\w.])(\w+)\s*\(", code):
        name = m.group(1)
        if name in KEYWORDS:
            continue
        line = code[: m.start()].count("\n") + 1
        calls[name].append(line)

    # ★★ 局部函数的**形参名**必须单独收集排除。
    #   理由：`M(r, :)` / `A(:)` 这类**索引**在正则看来就是"函数调用"，
    #   若形参名不在排除集合里，会被误报成"未定义/不可见函数"。
    #   （与 `@(A)` 匿名函数形参同一个坑。）
    params = set()
    for m in re.finditer(r"^[ \t]*function\b[^\n]*?\(([^)]*)\)", raw, re.M):
        for p in m.group(1).split(","):
            p = p.strip().lstrip("~").split("=")[0].strip()
            if re.fullmatch(r"[A-Za-z]\w*", p or ""):
                params.add(p)

    return {
        "path": path,
        "name": os.path.splitext(os.path.basename(path))[0],
        "code": code,
        "nocmt": nocmt,
        "raw": raw,
        "local": local_funs,
        "assigned": assigned,
        "params": params,
        "calls": calls,
    }


def _strip_strings(line):
    """把 MATLAB 字符串字面量替换成空格，避免把字符串里的文字误判成代码。

    ★ 专门为"续行块语法检查"服务：`'a...b'` 这种字符串里也可能出现 `...`，
      若不去掉字符串，会把它误当成续行符。
    """
    out = []
    i = 0
    quote = None
    n = len(line)
    while i < n:
        c = line[i]
        if quote:
            if c == "'" and i + 1 < n and line[i + 1] == "'":   # MATLAB 转义 ''
                out.append(" ")
                i += 2
                continue
            if c == quote:
                quote = None
            out.append(" ")
            i += 1
            continue
        if c in ("'", '"'):
            quote = c
            out.append(" ")
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def struct_fields(text, open_paren_index):
    """Return the set of field names inside struct( ... ) with balanced parens."""
    depth = 0
    i = open_paren_index
    n = len(text)
    while i < n:
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    body = text[open_paren_index + 1:i]
    return set(re.findall(r"'(\w+)'\s*,", body))


def strip_comments_and_strings(code):
    """Remove `%...` comments and '...' / "..." string literals.

    ★ Why this is necessary: block_check() tokenises the whole file. If a
    comment or a wrapped error() message happens to contain the words
    "if" / "for" / "end" (very easy in Chinese prose that enumerates
    "(1) ... (2) ..." or quotes a formula), those words get counted as
    block openers/terminators and the balance report goes wrong.
    This produced a bogus "块起始 3 / end 2" FAIL on
    crazyflie_slung_parameters.m while the file was in fact balanced.
    """
    out = []
    for line in code.split("\n"):
        # MATLAB block comments %{ ... %} are rare here; handle line comments.
        # A '%' starts a comment unless escaped ('%%' inside sprintf is still
        # a comment start in MATLAB source, which is what we want).
        # Strip strings first so that '%' inside a string is not treated as
        # a comment start.
        res = []
        i = 0
        n = len(line)
        while i < n:
            ch = line[i]
            if ch == "'":
                # transpose operator vs string literal: a quote is a string
                # opener when the previous non-space char is not an identifier
                # char, ')', ']', '}', or '.'
                prev = ""
                for j in range(len(res) - 1, -1, -1):
                    if not res[j].isspace():
                        prev = res[j]
                        break
                if prev and (prev.isalnum() or prev in "_)]}."):
                    res.append(ch)          # transpose
                    i += 1
                    continue
                # string literal: consume to closing quote ('' = escaped quote)
                i += 1
                while i < n:
                    if line[i] == "'":
                        if i + 1 < n and line[i + 1] == "'":
                            i += 2
                            continue
                        i += 1
                        break
                    i += 1
                continue
            if ch == "%":
                break                       # rest of line is a comment
            res.append(ch)
            i += 1
        out.append("".join(res))
    return "\n".join(out)


def block_check(code):
    """Count block openers vs `end` at statement level, and bracket balance.

    `end` used as an index (x(end)) lives inside brackets, so restricting the
    count to bracket-depth 0 correctly separates it from block-terminating `end`.

    ★ `switch` is followed by its expression and then `case`/`otherwise`, so a
    `switch` must be paired with exactly one terminating `end` -- handled by the
    generic opener rule below. `case`/`otherwise` are NOT openers and must not
    add depth (they are branches of the pending switch).
    """
    depth = 0
    n_end = 0
    n_open = 0
    openers = ("function", "if", "for", "while", "switch", "try", "parfor", "spmd")
    for m in re.finditer(r"[()\[\]{}]|\b[A-Za-z_]\w*\b", code):
        t = m.group(0)
        if t in "([{":
            depth += 1
        elif t in ")]}":
            depth -= 1
        elif depth == 0:
            if t == "end":
                n_end += 1
            elif t in openers:
                n_open += 1
    return n_open, n_end, depth


def main():
    infos = [parse_file(os.path.join(V3, f)) for f in FILES]
    # ★ 同目录的其它 .m 文件也是可见的（MATLAB 按文件名查找函数），
    #   把它们登记进 byname，否则 parameters.m 调用 crazyflie_slung_reference()
    #   会被误报为"未定义/不可见"。
    sibling_names = {os.path.splitext(f)[0] for f in SIBLING_M_FILES}
    byname = {i["name"]: i for i in infos}
    byname.update({nm: {"path": os.path.join(V3, nm + ".m")} for nm in sibling_names})

    print("=" * 78)
    print(" 5) 变量使用顺序检查（先用后定义 —— MATLAB 运行期才报错的坑）")
    print("=" * 78)
    order_problems = 0
    for info in infos:
        found = use_before_assign(info["raw"], info["local"])
        # 过滤：只保留确实不在同作用域更早处定义的
        if found:
            print(f"\n  [{info['name']}.m]  疑似先用后定义:")
            for name, line in found:
                print(f"      - {name}  行 {line}")
                order_problems += 1
        else:
            print(f"  [OK] {info['name']}.m")
    print(f"\n  --> 使用顺序问题: {order_problems}")

    # ------------------------------------------------- nargout/nargin 分派
    print()
    print("=" * 78)
    print(" 6) nargout/nargin 分派检查（匿名句柄恒返回 -1 的坑）")
    print("=" * 78)
    dispatch_problems = 0
    for info in infos:
        risks = nargout_dispatch_risk(info["nocmt"])
        if risks:
            print(f"\n  [{info['name']}.m]  分派未处理 -1（匿名句柄）:")
            for fn, arg, var, line in risks:
                print(f"      - {fn}({arg}) 行 {line} -> 变量 {var} 无 '-1' 守卫")
                dispatch_problems += 1
        else:
            print(f"  [OK] {info['name']}.m")
    print(f"\n  --> 分派问题: {dispatch_problems}")

    print("=" * 78)
    print(" 1) 函数可见性检查（MATLAB 局部函数是文件私有的）")
    print("=" * 78)
    problems = 0
    for info in infos:
        unknown = []
        for name, lines in sorted(info["calls"].items()):
            if name in info["local"]:
                continue
            if name in byname:                       # another .m file
                continue
            if name in BUILTINS:
                continue
            if name in info["assigned"]:             # indexing a variable
                continue
            if name in info["params"]:               # 局部函数形参（M(r,:) 等索引会被误判）
                continue
            if name in CALLABLE_VARS or name.endswith("Fcn"):
                continue
            unknown.append((name, lines))
        if unknown:
            print(f"\n  [{info['name']}.m]  未定义/不可见:")
            for name, lines in unknown:
                print(f"      - {name}()  行 {lines}")
                problems += 1
        else:
            print(f"  [OK] {info['name']}.m")
    print(f"\n  --> 可见性问题: {problems}")

    # ---------------------------------------------------------------- cfg paths
    print()
    print("=" * 78)
    print(" 2) cfg.* 字段路径检查")
    print("=" * 78)
    p = byname["crazyflie_slung_parameters"]
    pcode = p["nocmt"]
    groups = defaultdict(set)
    for m in re.finditer(r"cfg\.(\w+)\s*=\s*struct\s*\(", pcode):
        groups[m.group(1)] |= struct_fields(pcode, m.end() - 1)
    for m in re.finditer(r"cfg\.(\w+)\.(\w+)\s*=", pcode):
        groups[m.group(1)].add(m.group(2))
    alltop = set(re.findall(r"cfg\.(\w+)\s*=(?!=)", pcode))
    scalars = alltop - set(groups)

    print("  参数文件中定义的组与字段:")
    for g in sorted(groups):
        print(f"    cfg.{g:18s} : {len(groups[g]):2d} 字段  {{{', '.join(sorted(groups[g]))}}}")
    print(f"  标量字段: {{{', '.join(sorted(scalars))}}}")

    bad = 0
    for info in infos:
        for m in re.finditer(r"cfg\.(\w+)\.(\w+)", info["code"]):
            g, f = m.group(1), m.group(2)
            if g in groups and f not in groups[g]:
                line = info["code"][: m.start()].count("\n") + 1
                print(f"  [缺失] {info['name']}.m:{line}  cfg.{g}.{f}")
                bad += 1
    for info in infos:
        for m in re.finditer(r"cfg\.(\w+)(?!\s*\.)(?!\w)", info["code"]):
            nm = m.group(1)
            if nm in groups or nm in scalars or nm in BUILTINS or nm in KEYWORDS:
                continue
            line = info["code"][: m.start()].count("\n") + 1
            print(f"  [可疑] {info['name']}.m:{line}  cfg.{nm}")
            bad += 1
    print(f"\n  --> cfg 路径问题: {bad}")

    # ------------------------------------------------- state/command field names
    print()
    print("=" * 78)
    print(" 3) state.* / command.* / memory.* / derivative.* / summary.* 字段一致性")
    print("=" * 78)
    registry = {k: set() for k in
                ("state", "command", "memory", "desired", "derivative", "summary", "sim")}
    for info in infos:
        t = info["nocmt"]
        for var in registry:
            for m in re.finditer(rf"\b{var}\.(\w+)\s*=", t):
                registry[var].add(m.group(1))
            for m in re.finditer(rf"\b{var}\s*=\s*struct\s*\(", t):
                registry[var] |= struct_fields(t, m.end() - 1)
    for var in sorted(registry):
        if registry[var]:
            print(f"    {var}.{{{', '.join(sorted(registry[var]))}}}")

    bad2 = 0
    for info in infos:
        code = info["code"]
        for var in registry:
            if not registry[var]:
                continue
            for m in re.finditer(rf"(?<![\w.]){var}\.(\w+)", code):
                f = m.group(1)
                if f not in registry[var]:
                    line = code[: m.start()].count("\n") + 1
                    print(f"  [未定义] {info['name']}.m:{line}  {var}.{f}")
                    bad2 += 1
    print(f"\n  --> 字段问题: {bad2}")

    print()
    print("=" * 78)
    print(" 4) 代码块 / 括号配平（缺失 end、缺失右括号）")
    print("=" * 78)
    bad3 = 0
    for info in infos:
        # ★ 必须先剥掉注释与字符串字面量：中文注释里出现的 "if"/"for"/"end"
        #   会被误当成块起始/结束，造成虚假 FAIL（parameters.m 曾因此误报）。
        t = strip_comments_and_strings(info["code"])
        n_open, n_end, depth = block_check(t)
        ok = (n_open == n_end) and (depth == 0)
        flag = "OK  " if ok else "FAIL"
        if not ok:
            bad3 += 1
        print(f"  [{flag}] {info['name']}.m : 块起始 {n_open:3d} / end {n_end:3d} "
              f"| 括号净深度 {depth}")
        if not t.rstrip().endswith("end") and not t.rstrip().endswith("..."):
            print(f"         ^ 文件未以 end 结尾，请确认最后一行是否完整")
    print(f"\n  --> 配平问题: {bad3}")

    # ------------------------------------------------- 7) 续行块语法检查
    # ★★★ 这一类的由来：曾把多行注释**插进** `struct(...)` 的 `...` 续行块中间，
    #     而续行块里**每一行（包括注释行）都必须以 `...` 结尾**。
    #     结果语句被提前截断，MATLAB 直接报
    #         "文件: crazyflie_slung_parameters.m 行: 549 列: 48 无效表达式"
    #     —— 这是**硬语法错误**，却逃过了前 6 类检查（前 6 类都不做续行语法分析）。
    #     本检查专门堵这个洞。
    print()
    print(" 7) 续行块语法检查（`...` 续行块内的注释行必须以 `...` 结尾）")
    print("=" * 78)
    bad4 = 0
    for info in infos:
        hits = []
        lines = info["raw"].split("\n")
        prev_cont = False
        for n0, l in enumerate(lines, 1):
            s = l.strip()
            cont = _strip_strings(l).rstrip().endswith("...")
            if prev_cont and s.startswith("%") and not cont:
                hits.append((n0, s[:60]))
            prev_cont = cont
        if hits:
            print(f"\n  [{info['name']}.m]  续行块被注释行截断（会导致『无效表达式』）:")
            for n0, s in hits:
                print(f"      - 行 {n0}: {s}")
                bad4 += 1
        else:
            print(f"  [OK] {info['name']}.m")
    print(f"\n  --> 续行语法问题: {bad4}")

    print()
    print("=" * 78)
    print(f" 汇总: 可见性 {problems} | 使用顺序 {order_problems} "
          f"| 分派 {dispatch_problems} | cfg路径 {bad} "
          f"| 字段 {bad2} | 配平 {bad3} | 续行 {bad4}")
    print("=" * 78)
    return 0 if (problems == 0 and order_problems == 0 and dispatch_problems == 0
                 and bad == 0 and bad2 == 0 and bad3 == 0 and bad4 == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
