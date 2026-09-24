function [derivative, info] = crazyflie_slung_dynamics(state, uAll, cfg)
%CRAZYFLIE_SLUNG_DYNAMICS 多架四旋翼协同吊运刚体负载的完整动力学。
%
% 实现 Lee 2014 (arXiv:1403.3684) / Lee 2018 (IEEE TCST) 的式 (5)-(8)，
% 对应 2014 版的式 (6)-(9)：
%
%   (5) 负载平动（n 机耦合）
%       Mq (x0_ddot - g e3) - sum_i m_i q_i q_i' R0 hat(rho_i) Omega0_dot
%           = sum_i u||_i - sum_i m_i l_i ||omega_i||^2 q_i - sum_i m_i q_i q_i' R0 hat(Omega0)^2 rho_i
%   (6) 负载转动
%       ( J0 - sum_i m_i hat(rho_i) R0' q_i q_i' R0 hat(rho_i) ) Omega0_dot
%         + sum_i m_i hat(rho_i) R0' q_i q_i' (x0_ddot - g e3) + Omega0 x J0 Omega0
%           = sum_i hat(rho_i) R0' ( u||_i - m_i l_i ||omega_i||^2 q_i - m_i q_i q_i' R0 hat(Omega0)^2 rho_i )
%   (7) 第 i 根绳索
%       omega_dot_i = (1/l_i) hat(q_i) a_i - (1/(m_i l_i)) hat(q_i) u_perp_i
%   (8) 第 i 架四旋翼姿态
%       J_i Omega_dot_i + Omega_i x J_i Omega_i = M_i
%
% 其中等效质量矩阵（论文 (5) 之后正文）
%       Mq = m0 I + sum_i m_i q_i q_i'
% ★ 注意 Mq 的对角系数是 **m0**（负载质量），不是 (m0 + sum m_i)。
%   一致性校验：把 (17) u||_i = mu_i + m_i l_i||omega_i||^2 q_i + m_i q_i q_i' a_i
%   代入 (5)，耦合项相消后可得论文 (18)：m0(x0_ddot - g e3) = sum_i mu_i，
%   即绳索张力之和恰好平衡负载自身重力。若把 Mq 的对角写成 (m0+sum m_i)，
%   这一恒等式就不成立，且无法与 (18) 对应。
%
% ★★ 关于 a_i = x0_ddot - g e3 - alpha_i（alpha_i := R0 hat(rho_i) Omega0_dot）
%   论文把 a_i 当作已知量写进 (16)(17)(27)，但 (5)(6) 表明 x0_ddot 与 Omega0_dot
%   本身又依赖 u||_i，而 u||_i 又依赖 a_i —— 这是一个代数环。论文的处理方式是在
%   推导 (5)(6) 时把 (17) 代入，从而消去该环。本文件沿用同一做法：先把 (5)(6)
%   整理成关于 [X; Omega0_dot] 的 6x6 线性系统解出二者，再回代算出 a_i。
%
%   ⚠ 曾经尝试把 alpha_i 当作**独立未知量**再补一组定义方程 alpha_i = B_i Omega0_dot
%     来解同一个环。这是错的：那样 (5) 会通过 alpha 拿到一份 Omega0_dot 的贡献，
%     而 (6) 又直接含 Omega0_dot，同一自由度被重复计入。判别方法：把解代回论文 (18)
%        m0 (x0_ddot - g e3) = sum_i mu_i
%     6x6 消元形式残差 2.3e-5 N；alpha 独立形式残差 2.4e-2 N（差三个数量级，物理错误）。
%
% 状态结构体字段：
%   loadPosition (3)         x0
%   loadVelocity (3)         x0_dot
%   loadRotation (3x3)       R0
%   loadBodyRate (3)         Omega0
%   linkUnits (3xn)          第 i 列 = q_i（由第 i 架四旋翼指向负载）
%   linkRates (3xn)          第 i 列 = q_i_dot
%   rotations (3x3xn)        第 i 页 = R_i
%   bodyRates (3xn)          第 i 列 = Omega_i
%   bodyTorques (3xn)        第 i 列 = M_i（由仿真主程序的内环速率环给出）
%
% 输出额外字段：
%   linkPerpForces (3xn)     第 i 列 = u_perp_i，供仿真主程序的诊断项
%                            l*||hat(q_i) u_perp_i|| 使用（避免因
%                            u_perp_i 含 u|| 分量而把它记为"结构误差"）

% ------------------------------------------------------------------ 状态读取
x0 = state.loadPosition(:);
v0 = state.loadVelocity(:);
R0 = projectSO3(state.loadRotation);
Omega0 = state.loadBodyRate(:);
Qall = state.linkUnits;                       % 3 x n
QDall = state.linkRates;                      % 3 x n
Rall = state.rotations;                       % 3 x 3 x n
OmAll = state.bodyRates;                      % 3 x n

% ------------------------------------------------------------------ 参数读取
e3 = [0; 0; 1];
g = cfg.vehicle.gravity;
m0 = cfg.payload.mass;
J0 = cfg.payload.inertia;
rhoAll = cfg.payload.attachPoints;            % 3 x n
m = cfg.vehicle.mass;
J = cfg.vehicle.inertia;
l = cfg.link.length;
n = cfg.vehicle.count;
uAll = reshape(uAll, 3, n);

% hat(Omega0)^2 * rho_i 与 hat(rho_i)（对所有 i 通用）
hatOm0 = hat(Omega0);
hatOm0Sq = hatOm0 * hatOm0;

% ---------------------------------------------- 逐机累加，装配耦合系统
% 未知量 Z_p = [ X ; Omega0_dot ]，X = x0_ddot - g e3（论文写法）
sumQ = zeros(3, 3);          % sum_i m_i q_i q_i'
sumMRho = zeros(3, 3);       % sum_i m_i q_i q_i' R0 hat(rho_i)          -> (1,2) 块（取负）
sumRhoQ = zeros(3, 3);       % sum_i m_i hat(rho_i) R0' q_i q_i'        -> (2,1) 块
sumRhoQRho = zeros(3, 3);    % sum_i m_i hat(rho_i) R0' q_i q_i' R0 hat(rho_i)
rhs1 = zeros(3, 1);          % (5) 右端（不含 Omega0_dot 耦合项）
rhs2 = zeros(3, 1);          % (6) 右端（不含 -Omega0 x J0 Omega0）

uParallelAll = zeros(3, n);
uPerpAll = zeros(3, n);
omegaAll = zeros(3, n);

for i = 1:n
    qi = normalizeVector(Qall(:, i));
    qdi = QDall(:, i);
    rhoi = rhoAll(:, i);

    % 绳索角速度（论文 (1)：q_dot = omega x q，故 omega = q x q_dot）
    omegai = cross(qi, qdi);
    omegaAll(:, i) = omegai;

    % 控制力按论文 (9)(10) 分解：u|| = q q' u，u_perp = (I - q q') u
    ui = uAll(:, i);
    uPar = dot(qi, ui) * qi;
    uPerp = ui - uPar;
    uParallelAll(:, i) = uPar;
    uPerpAll(:, i) = uPerp;

    Qi = qi * qi.';
    rhoHat = hat(rhoi);
    hatOm0SqRho = hatOm0Sq * rhoi;

    % (5) 中第 i 机的贡献
    zeta_i = uPar - m * l * dot(omegai, omegai) * qi - m * Qi * R0 * hatOm0SqRho;
    rhs1 = rhs1 + zeta_i;

    % (6) 中第 i 机的贡献
    % 论文 (6) 的力矩臂项为 hat(rho_i) R0' zeta_i。
    % 这里不能写成 rhoHat.'：那会把分配矩阵 P 中的期望力矩反号，
    % 尤其在 yaw 通道开启时形成正反馈。
    rhs2 = rhs2 + rhoHat * R0.' * zeta_i;

    % 各耦合矩阵块
    sumQ = sumQ + m * Qi;
    sumMRho = sumMRho + m * Qi * R0 * rhoHat;
    R0tQi = R0.' * Qi;
    sumRhoQ = sumRhoQ + m * rhoHat * R0tQi;
    sumRhoQRho = sumRhoQRho + m * rhoHat * R0tQi * R0 * rhoHat;
end

% 论文 (5)：Mq = m0 I + sum_i m_i q_i q_i'
Mq = m0 * eye(3) + sumQ;

% ----------------------------------------------------------------------
% 联立 (5)(6)，未知量 Z_p = [ X ; Omega0_dot ]
%   [ Mq               , -sumMRho          ] [ X          ]   [ rhs1 ]
%   [ sumRhoQ          ,  J0 - sumRhoQRho  ] [ Omega0_dot ] = [ rhs2 - Omega0xJ0Omega0 ]
%
% ★ 这里必须保持 6x6 的**消元后**形式：Omega0_dot 是唯一未知量，其系数矩阵
%   由 R0 hat(rho_i) 这样的**常量**构成。曾尝试把 alpha_i := R0 hat(rho_i) Omega0_dot
%   作为独立未知量再补一组定义方程 alpha_i = R0 hat(rho_i) Omega0_dot，
%   结果 (5) 与 (6) 会各自得到一份 Omega0_dot 的贡献（(5) 通过 alpha、(6) 直接），
%   同一自由度被重复计入，解出的 x0_ddot 不再满足论文 (18)
%   校验方式：把解代回 m0(x0_ddot - g e3) = sum_i mu_i，
%   6x6 形式残差 2.3e-5 N，alpha 独立形式残差 2.4e-2 N（差三个数量级）。
% ----------------------------------------------------------------------
Mblock = [Mq,        -sumMRho; ...
          sumRhoQ,   J0 - sumRhoQRho];
% ★★★ 负载**转动阻尼**：作为**外部力矩**加在转动方程右端（不是控制指令）。
%   转成控制需求会经 "Md -> 张力分配 -> q_id 水平化" 形成正反馈（实测必发散），
%   而外部力矩直接作用在负载上，绕开控制器 ⇒ 稳定。
%   物理来源与取值依据见 crazyflie_slung_parameters.m 中 rotationalDamping 的说明。
%   作用：压掉"负载偏航 <-> 绳索扭转"那个既不可控又无阻尼的模态。
dampingTorque = -cfg.payload.rotationalDamping * Omega0;
rhsBlock = [rhs1; rhs2 - cross(Omega0, J0 * Omega0) + dampingTorque];

% ★★ 发散保护 + **一次性详细报告**
%   ⚠ 下面这个分支是"让故障可读"，**不是把故障修好**。
%   它在矩阵已含 NaN/Inf 时把解置零，于是状态不会继续被 NaN 污染 ——
%   但这也意味着**状态日志可能全部保持有限**，用"日志里有没有 NaN"查不出问题。
%   因此这里必须把**病因**直接打出来。
%   ★ 判据用 `isnan`/`isinf`（base MATLAB 一定有），不用 `isfinite`。
%     （也不写成 `@(A)` 匿名函数：静态检查会把匿名参数 A 误认成函数调用。）
%   ★★ 踩过的坑：`warning` 必须放在"只报一次"分支**里面**。
%      放在外面会变成"每个坏步都刷一条" —— 坏步有 2 万步时，
%      真正的病因报告会被彻底冲掉（实测刷屏 20812 条）。
nBadM = nnz(isnan(Mblock)) + nnz(isinf(Mblock));
nBadR = nnz(isnan(rhsBlock)) + nnz(isinf(rhsBlock));
finiteInputs = (nBadM == 0) && (nBadR == 0);
if finiteInputs
    solution = Mblock \ rhsBlock;
else
    persistent nDynCall reportedNonFinite
    if isempty(nDynCall)
        nDynCall = 0;
        reportedNonFinite = false;
    end
    nDynCall = nDynCall + 1;
    if ~reportedNonFinite
        reportedNonFinite = true;
        % 报告同时写到**文件**与**命令行**：
        %   坏步可能上万条，命令行输出容易被冲掉；但用户又需要能直接复制。
        %   做法：先写文件，再把文件内容回放到命令行。
        reportPath = fullfile(tempdir, 'cf_slung_nonfinite_report.txt');
        fid = fopen(reportPath, 'w');
        if fid < 0
            fid = 1;        % 打不开就只写命令行
        end
        fprintf(fid, '*** 6x6 代数系统出现非有限量（这是**第 %d 个坏步**，不是第 %d 次调用）***\n', ...
            nDynCall, nDynCall);
        fprintf(fid, '    ★ 判据只检查 Mblock/rhsBlock 本身；若 Mblock 全 ok 而 rhsBlock 全 NaN，\n');
        fprintf(fid, '      说明 NaN 来自**上游**（控制器/上一拍的解），本报告下面的 dump 用于定位那个上游量。\n');
        fprintf(fid, '    矩阵规模 Mblock = %dx%d (%s), rhsBlock = %dx%d (%s)\n', ...
            size(Mblock, 1), size(Mblock, 2), class(Mblock), ...
            size(rhsBlock, 1), size(rhsBlock, 2), class(rhsBlock));
        fprintf(fid, '    非有限计数：Mblock 中共 %d 个（NaN %d / Inf %d）\n', ...
            nBadM, nnz(isnan(Mblock)), nnz(isinf(Mblock)));
        fprintf(fid, '                rhsBlock 中共 %d 个（NaN %d / Inf %d）\n', ...
            nBadR, nnz(isnan(rhsBlock)), nnz(isinf(rhsBlock)));
        fprintf(fid, '    逐项报告（**第一行标 BAD 的就是病因**；全 ok 则看下面的原始 dump）\n');
        reportFinite(fid, 'rhs1', rhs1);
        reportFinite(fid, 'rhs2', rhs2);
        reportFinite(fid, 'Mq', Mq);
        reportFinite(fid, 'sumMRho', sumMRho);
        reportFinite(fid, 'sumRhoQ', sumRhoQ);
        reportFinite(fid, 'sumRhoQRho', sumRhoQRho);
        reportFinite(fid, 'Omega0 (负载角速度)', Omega0);
        reportFinite(fid, 'R0 (负载姿态)', R0);
        reportFinite(fid, 'Qall (绳向 q_i)', Qall);
        reportFinite(fid, 'QDall (绳向速率)', QDall);
        reportFinite(fid, 'omegaAll (绳向角速度)', omegaAll);
        reportFinite(fid, 'uAll (实际作用力 ★唯一非状态输入)', uAll);
        reportFinite(fid, 'v0 (负载速度)', v0);
        reportFinite(fid, 'hatOm0Sq', hatOm0Sq);
        reportFinite(fid, 'rhoAll (挂点几何)', rhoAll);
        reportFinite(fid, '标量 m0/m/l/g', [m0; m; l; g]);
        reportFinite(fid, 'J0 (负载惯量)', J0);
        fprintf(fid, '    ---- Mblock (6x6) ----\n');
        printMatrix(fid, Mblock);
        fprintf(fid, '    ---- rhsBlock (6x1) ----\n');
        printMatrix(fid, rhsBlock);
        fprintf(fid, '    ---- uAll (3x%d) ----\n', n);
        printMatrix(fid, uAll);
        fprintf(fid, '    ---- Qall (绳向 q_i, 3x%d) ----\n', n);
        printMatrix(fid, Qall);
        fprintf(fid, '    ---- QDall (绳向速率, 3x%d) ----\n', n);
        printMatrix(fid, QDall);
        fprintf(fid, '    ---- R0 (负载姿态) ----\n');
        printMatrix(fid, R0);
        fprintf(fid, '    ---- Omega0 (负载角速度) ----\n');
        printMatrix(fid, Omega0);
        fprintf(fid, '    --------------------------------------------------------\n');
        if fid ~= 1
            fclose(fid);
            % 回放到命令行，方便用户直接复制
            try
                fprintf('%s', fileread(reportPath));
                fprintf('    （同一份报告也已写入：%s）\n\n', reportPath);
            catch
                fprintf('    报告已写入：%s\n\n', reportPath);
            end
        end
        warning('crazyflie_slung:NonFiniteState', ...
            '检测到 6x6 代数系统出现非有限量；详细报告已一次性打印（本警告只出一次）。');
    end
    solution = zeros(6, 1);
end
% ★ 把"这一步是否是坏步"报出去，供仿真主程序累计成 summary.nonFiniteSolveCount。
%   这个计数是**唯一不会骗人的判据**：状态日志可能全有限
%   （因为坏步的解被置零、NaN 没传进状态），只看日志会误判"没有发散"。
info.nonFiniteInputs = ~finiteInputs;
X = solution(1:3);
x0DotDot = X + g * e3;
Omega0Dot = solution(4:6);

% ------------------------------------------------ 逐机绳索与姿态动力学
qDotDotAll = zeros(3, n);
omegaDotAll = zeros(3, n);
OmegaDotAll = zeros(3, n);

for i = 1:n
    qi = normalizeVector(Qall(:, i));
    omegai = omegaAll(:, i);
    rhoi = rhoAll(:, i);
    rhoHat = hat(rhoi);

    % 论文 (16)：挂点相对加速度
    %   a_i = x0_ddot - g e3 + R0 hat(Omega0)^2 rho_i - R0 hat(rho_i) Omega0_dot
    %      = X - alpha_i,  alpha_i := R0 hat(rho_i) Omega0_dot
    alphaI = R0 * rhoHat * Omega0Dot;
    ai = X - alphaI;

    % 论文 (7)：绳索角加速度
    omegaDotI = (1 / l) * cross(qi, ai) - (1 / (m * l)) * cross(qi, uPerpAll(:, i));
    omegaDotAll(:, i) = omegaDotI;

    % 由 q_ddot = omega_dot x q + omega x (omega x q) 得到绳索向量的二阶导
    qDotDotAll(:, i) = cross(omegaDotI, qi) + cross(omegai, cross(omegai, qi));

    % 论文 (8)：机体姿态动力学。
    % 本方案的角速度内环在仿真主程序中把 Omega_cmd_i 转成等效力矩 M_i，
    % 此处只负责由力矩推进机体角速度。
    OmegaI = OmAll(:, i);
    Mi = state.bodyTorques(:, i);
    OmegaDotAll(:, i) = J \ (Mi - cross(OmegaI, J * OmegaI));
end

% ------------------------------------------------------------------ 输出
derivative = struct();
derivative.loadPosition = v0;
derivative.loadVelocity = x0DotDot;
derivative.loadRotation = R0 * hat(Omega0);
derivative.loadBodyRate = Omega0Dot;
derivative.linkUnits = QDall;
derivative.linkRates = qDotDotAll;
derivative.rotations = zeros(3, 3, n);
for i = 1:n
    derivative.rotations(:, :, i) = Rall(:, :, i) * hat(OmAll(:, i));
end
derivative.bodyRates = OmegaDotAll;
% 中间量（供仿真主程序与日志使用）
derivative.loadAcceleration = x0DotDot;
derivative.loadBodyAcceleration = Omega0Dot;
derivative.linkAngularVelocities = omegaAll;
derivative.linkAngularAccelerations = omegaDotAll;
derivative.parallelForces = uParallelAll;
% ★ 这里输出**原始** u_perp_i（由 (I - q q') u_i 得到）。诊断项必须用
%   l * ||hat(q_i) u_perp_i||，因为 (7) 只用到 hat(q_i) u_perp_i，
%   u_perp_i 沿 q_i 的分量对绳索动力学无影响，不应被记为"结构误差"。
derivative.linkPerpForces = uPerpAll;
derivative.perpendicularForces = uPerpAll;
derivative.massMatrix = Mq;
end

% ======================================================================
function S = hat(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function q = normalizeVector(q)
q = q(:);
n = norm(q);
if n < eps
    q = [0; 0; 1];      % 退化时取 +e3
else
    q = q / n;
end
end

function R = projectSO3(R)
[U, ~, V] = svd(R);
R = U * V.';
if det(R) < 0
    U(:, 3) = -U(:, 3);
    R = U * V.';
end
end

% ======================================================================
function reportFinite(fid, name, M)
% 逐项报告"这个量里有没有非有限值"。用于发散时一次性定位病因。
% ★ 除了线性位置，还打出 **(行, 列)** —— 对 6x6 的 Mblock 来说，
%   行 1~3 = (5) 的平动方程、行 4~6 = (6) 的转动方程；
%   列 1~3 = X 的系数、列 4~6 = Omega0_dot 的系数。看到行列就能定位是哪一块坏。
% 只打印前 8 个越界位置，避免刷屏。
bad = isnan(M) | isinf(M);
nb = nnz(bad);
if nb == 0
    fprintf(fid, '    [ok  ] %s\n', name);
else
    idx = find(bad(:).');
    nShow = min(numel(idx), 8);
    [rr, cc] = ind2sub(size(M), idx(1:nShow));
    pairs = arrayfun(@(k) sprintf('(%d,%d)', rr(k), cc(k)), 1:nShow, ...
        'UniformOutput', false);
    if nb == numel(M)
        what = '全部非有限';
    elseif nnz(isnan(M)) == nb
        what = '全是 NaN';
    elseif nnz(isinf(M)) == nb
        what = '全是 Inf';
    else
        what = sprintf('NaN %d / Inf %d', nnz(isnan(M)), nb - nnz(isnan(M)));
    end
    fprintf(fid, '    [BAD ] %-28s %d/%d 非有限（%s）；位置 %s\n', ...
        name, nb, numel(M), what, strjoin(pairs, ' '));
end
end

% ======================================================================
function printMatrix(fid, M)
% 把矩阵完整打到 fid（NaN/Inf 照原样显示，便于肉眼比对）。
for r = 1:size(M, 1)
    fprintf(fid, '   ');
    fprintf(fid, '%12.5g ', M(r, :));
    fprintf(fid, '\n');
end
end
