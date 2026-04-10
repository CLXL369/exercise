function lam = find_lambda(r, VN0, VC, Q_c, G, L, hb2_2mu)
% FIND_LAMBDA  通过薛定谔本征值条件求归一化因子 lambda。
%
%   lam = find_lambda(r, VN0, VC, Q_c, G, L, hb2_2mu)
%
% =========================================================================
% 方法 — 薛定谔本征值条件（Bai & Ren 2018，第2节）
% =========================================================================
%
%   文章指出："归一化因子 lambda 由要求 alpha-核有效势在 L=0 时
%   重现实验测量的 Q_alpha 值来确定。"
%
%   即：求 lambda，使得径向薛定谔方程
%
%       [ -hbar^2/(2*mu) * d^2/dr^2  +  V_eff(r) ] u(r)  =  Q_c * u(r)
%
%   在内部经典允许区 [0, r1] 内恰好有 N_target = (G-L)/2 = 8 个节点，
%   且解是可归一化的。
%   这是用数值积分精确求解的 Wildermuth 量子化条件（非 WKB 近似）。
%
%   算法 — 节点计数二分法（Sturm-Liouville 定理）：
%     1. 固定能量 E = Q_c（实验中心值）。
%     2. 用 Numerov 方法从 r=0 向外积分薛定谔方程。
%     3. 计算 u(r) 在 [0, r1]（内部经典允许区）中的节点数。
%     4. 增大 lambda 使势阱加深 → 在能量 E=Q_c 处容纳更多节点。
%     5. 每当一个新节点进入势阱，u 在 r1 之前的符号翻转，
%        据此可精确二分确定 lambda*。
%
%   输入参数
%     r        – 径向网格（fm），列向量，从 dr 开始
%     VN0      – V_N(r)/lambda：未缩放核势（MeV）
%     VC       – V_C(r)：库仑势，与 lambda 无关（MeV）
%     Q_c      – 实验 Q_alpha 中心值（MeV）
%     G        – Wildermuth 全局量子数（= 16）
%     L        – 轨道角动量（基态跃迁取 0）
%     hb2_2mu  – hbar^2/(2*mu)，单位 MeV·fm^2
%   输出参数
%     lam – 满足 E_{N_target}(lam) = Q_c 的 lambda（薛定谔条件）

dr       = r(2) - r(1);
N_target = (G - L) / 2;    % 节点数 = 8（G=16, L=0 时，公式18）

% ---- 粗略扫描，定位区间 [lam_lo, lam_hi] ---------------------------------
% 对每个 lambda，统计 u(r; E=Q_c) 在内部经典允许区的节点数。
% 节点数从 N_target-1 变为 N_target 的位置即为包含 lambda* 的区间。
lam_arr  = linspace(0.25, 1.25, 60);
node_arr = zeros(1, numel(lam_arr));
u_arr    = zeros(1, numel(lam_arr));   % u 在势阱最后一点的值（仅符号有意义）

for i = 1:numel(lam_arr)
    V_eff = lam_arr(i) * VN0 + VC;
    if L > 0
        V_eff = V_eff + hb2_2mu * L*(L+1) ./ r.^2;
    end
    [node_arr(i), u_arr(i)] = nodes_and_u(r, V_eff, Q_c, hb2_2mu, dr);
end

% 扫描中节点数首次达到 N_target 的索引
idx_hi = find(node_arr >= N_target, 1);
if isempty(idx_hi)
    lam = lam_arr(end);
    warning('find_lambda: N_target=%d 个节点未达到；返回 lam=%.4f', N_target, lam);
    return;
end
if idx_hi == 1
    lam = lam_arr(1);
    warning('find_lambda: 在 lam=%.3f 以下未找到区间；返回 lam=%.4f', lam_arr(1), lam);
    return;
end

lam_lo = lam_arr(idx_hi - 1);
lam_hi = lam_arr(idx_hi);
u_lo   = u_arr(idx_hi - 1);

% ---- 对势阱最后一点 u 的符号做二分法 ------------------------------------
% lam_lo 处：N_target-1 个节点 → 符号为 u_lo
% lam_hi 处：N_target   个节点 → 符号与 u_lo 相反
% 二分直到 |lam_hi - lam_lo| < 1e-8
for iter = 1:80
    lam   = 0.5 * (lam_lo + lam_hi);
    V_eff = lam * VN0 + VC;
    if L > 0
        V_eff = V_eff + hb2_2mu * L*(L+1) ./ r.^2;
    end
    [~, u_m] = nodes_and_u(r, V_eff, Q_c, hb2_2mu, dr);

    if u_m * u_lo >= 0
        lam_lo = lam;
        u_lo   = u_m;
    else
        lam_hi = lam;
    end

    if (lam_hi - lam_lo) < 1e-8
        break;
    end
end

lam = 0.5 * (lam_lo + lam_hi);
end

% =========================================================================
function [n, u_last] = nodes_and_u(r, V_eff, E, hb2_2mu, dr)
% 计算 Numerov 解 u(r; E) 在经典允许内区 [0, r1] 中的节点数，
% 并返回 r1 之前最后一个网格点的 u 值。
%
% 化简径向薛定谔方程：
%   u''(r) + g(r)*u(r) = 0,   g(r) = (E - V_eff(r)) / hb2_2mu
%   u(0) = 0,  u(dr) = dr   （L=0 时原点正则边界条件）
%
% Numerov 递推格式（4阶精度，步长 h = dr）：
%   u_{n+1} = [ 2*u_n*(1 - 5h^2/12 * g_n) - u_{n-1}*(1 + h^2/12 * g_{n-1}) ]
%             / (1 + h^2/12 * g_{n+1})

N = length(r);
g = (E - V_eff) / hb2_2mu;    % g(r)：在势阱内（经典允许区）为正

% 内转折点：g 首次变为非正的索引
idx_r1 = find(g <= 0, 1);
if isempty(idx_r1) || idx_r1 <= 2
    n = 0;  u_last = r(1);  return;
end
n_steps = idx_r1 - 1;          % 从 r(1) 积分到 r(n_steps)

% r=0 处的 g 值（折叠势在原点处 V_eff(0) ≈ 0）
g0     = E / hb2_2mu;
u_prev = 0;                     % u(r = 0)
u_curr = r(1);                  % u(r = dr) = dr  [L=0 正则初始条件]
n      = 0;
s_prev = +1;                    % sign(u_curr) > 0

for j = 1 : n_steps - 1
    % 上一步、当前步、下一步的 g 值
    % 注意：j=1 时不能用 g(j-1)=g(0)（MATLAB 下标从1开始），
    %       改用 r=0 处的解析值 g0。
    if j == 1
        g_prev = g0;
    else
        g_prev = g(j - 1);
    end
    g_curr = g(j);
    g_next = g(j + 1);

    denom  = 1.0 + (dr*dr / 12.0) * g_next;
    u_next = (2.0 * u_curr * (1.0 - (5.0/12.0)*dr*dr*g_curr) ...
              - u_prev * (1.0 + (dr*dr/12.0)*g_prev)) / denom;

    % 统计符号变化次数（= 节点数）
    s_next = sign(u_next);
    if s_next ~= 0 && s_next * s_prev < 0
        n      = n + 1;
        s_prev = s_next;
    end

    % 每步重新归一化，防止浮点溢出（符号保留，仅幅度缩放）
    a      = abs(u_next);
    if a > 1e10
        scale  = a;
        u_prev = u_curr / scale;
        u_curr = u_next / scale;
    else
        u_prev = u_curr;
        u_curr = u_next;
    end
end

u_last = u_curr;    % 内转折点之前的 u 值（符号）
end
