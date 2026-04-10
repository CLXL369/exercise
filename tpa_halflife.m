function T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
% TPA_HALFLIFE  用双势方法（TPA）计算 alpha 衰变半衰期。
%
%   T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
%
% =========================================================================
% 理论 — 双势方法（TPA），Bai & Ren (2018) 第2节
% =========================================================================
%
% 在道半径 R（V 的外侧极大值处）将 alpha-核总势分为两个区域：
%
%   区域I  （内侧，r <= R）：V1(r) = V_N(r) + V_C(r)   [完整势]
%   区域II （外侧，r > R ）：V2(r) = V_C(r)             [仅库仑势]
%
% ---- 公式(7)：衰变宽度 -----------------------------------------------
%   Gamma_alpha = (4 * hbar^2 * ktilde^2) / (mu * k)
%                 * |phi_L(R)|^2
%                 * |chi_L(eta, kR)|^2
%
%   其中（所有量均采用 fm、MeV 单位）：
%     k      = sqrt(2*mu*Q_alpha) / hbar          渐近波数
%     ktilde = sqrt(2*mu*(V(R)-Q_alpha)) / hbar   R 处的虚波数
%     phi_L(R) = R 处的内区波函数，归一化条件 integral_0^R |phi|^2 dr = 1
%     chi_L  = F_L(eta, kR) = 正则库仑波函数 [公式(6) 区域II]
%     eta    = Z_alpha * Z_core * e^2 * mu / (hbar^2 * k)   [Sommerfeld 参数]
%
%   半衰期：T_{1/2} = hbar * ln(2) / Gamma_alpha        [P_alpha = 1]
%
% ---- |phi_L(R)|^2：R 处的 WKB 内区波函数 -------------------------
%   不使用从 r=0 直接向外的 Numerov 积分，因为它会捕获内势垒 [r1, R]
%   中指数增长的 WKB 模式，使 |phi(R)| 高估 exp(G_inner) 倍。
%
%   标准 TPA 实现改用 WKB 连接公式，通过内势垒匹配内势阱解：
%
%       |phi_L(R)|^2 = exp(-2 * G_inner) / (2 * ktilde * I_inv)
%
%   其中：
%     r1      = 内转折点：V(r) = Q_alpha 的第一个 r
%     G_inner = integral_{r1}^{R}  kappa(r)  dr    [内势垒作用量]
%     kappa(r)= sqrt((V(r)-Q_alpha)/hb2_2mu)        [势垒中局部衰减率]
%     I_inv   = integral_0^{r1}  1/k(r)  dr         [WKB 归一化积分]
%     k(r)    = sqrt((Q_alpha-V(r))/hb2_2mu)         [内势阱中的波数]
%
%   由于这些系统的 Sommerfeld 参数 eta >> 1，WKB 近似精度很高，
%   与公式(7) 用 WKB 方法计算 phi_L 等价。
%
% =========================================================================
%   所有长度单位：fm；能量单位：MeV。

dr = r(2) - r(1);

%% ---- 步骤1：定位道半径 R --------------------------------
% R 是 V(r) 在 r > 5 fm 之后的最外侧极大值处（公式7 中的道半径）。
% 通过检测 dV/dr 在 5 fm 之外最后一次从正变负来确定。
dV       = diff(V) / dr;
r_m      = 0.5*(r(1:end-1) + r(2:end));
mask     = r_m > 5.0;
dV_m     = dV(mask);
idx_mask = find(mask);

sgn_chg = find(dV_m(1:end-1) > 0 & dV_m(2:end) < 0, 1, 'last');
if ~isempty(sgn_chg)
    idx_R = idx_mask(sgn_chg);
else
    sub_idx = find(r > 5.0);
    [~, loc] = max(V(sub_idx));
    idx_R    = sub_idx(loc);
end
R = r(idx_R);

if V(idx_R) <= Q_a
    warning('tpa_halflife: V(R) <= Q_a；未找到势垒顶部。');
    T_half = Inf;
    return;
end

%% ---- 步骤2：R 处的波数（公式7）---------------------------------
% k      = sqrt(2*mu*Q_a)          / hbar    [渐近，区域II]
% ktilde = sqrt(2*mu*(V(R)-Q_a))   / hbar    [虚数，势垒顶部]
k      = sqrt(2 * mu_c2 * Q_a)              / hbar_c;
ktilde = sqrt(2 * mu_c2 * (V(idx_R) - Q_a)) / hbar_c;

%% ---- 步骤3：内转折点 r1 --------------------------
% r1 是 V(r) 首次超过 Q_alpha 的径向点（从原点向外），
% 即 [0, R] 区间内 V(r) = Q_alpha 的第一个根。
idx_r1 = find(V(1:idx_R) > Q_a, 1);
if isempty(idx_r1) || idx_r1 <= 1
    warning('tpa_halflife: 未找到内转折点。');
    T_half = Inf;
    return;
end

%% ---- 步骤4：WKB 归一化积分 I_inv = int_0^{r1} k^{-1} dr --
% k(r) = sqrt((Q_a - V(r)) / hb2_2mu)  在内势阱经典允许区。
% 被积函数 1/k(r) 有界（这些系统中 k >= ~0.25 fm^-1），
% 因此简单矩形求和即可给出可靠结果，无需奇点处理。
k2_in  = (Q_a - V(1:idx_r1-1)) / hb2_2mu;   % 在势阱内为正
k_in   = sqrt(max(k2_in, 0));
k_safe = max(k_in, 1e-6);                     % 安全下限（典型 Q 值下不会触发）
I_inv  = sum(1.0 ./ k_safe) * dr;             % fm^2（单位：[fm^-1]^-1 * fm = fm^2 ✓）

if I_inv <= 0
    warning('tpa_halflife: I_inv <= 0；无法归一化内区波函数。');
    T_half = Inf;
    return;
end

%% ---- 步骤5：内势垒作用量 G_inner = int_{r1}^{R} kappa dr -----
% kappa(r) = sqrt((V(r)-Q_a)/hb2_2mu)  在经典禁戒区。
kap_sq  = (V(idx_r1:idx_R) - Q_a) / hb2_2mu;
kap     = sqrt(max(kap_sq, 0));
G_inner = sum(kap) * dr;

%% ---- 步骤6：R 处的 WKB 内区波函数幅度 --------------------
% |phi_L(R)|^2 = exp(-2*G_inner) / (2 * ktilde * I_inv)   [恒为正]
uR_sq = exp(-2 * G_inner) / (2 * ktilde * I_inv);

%% ---- 步骤7：Sommerfeld 参数与库仑函数（公式6）--------
% eta = Z_a*Z_c*e^2*mu / (hbar^2*k)  [无量纲 Sommerfeld 参数]
% chi = F_L(eta, kR)                 [正则库仑函数，区域II]
eta = Z_a * Z_c * e2 * mu_c2 / (hbar_c^2 * k);
chi = coulomb_F(L, eta, k * R);

%% ---- 步骤8：衰变宽度与半衰期（公式7）-------------------------
% Gamma_alpha [MeV] = 4*(hbar^2*ktilde^2/(mu*k)) * |phi(R)|^2 * |chi|^2
%   单位验证：MeV^2·fm^2 * fm^-2 / (MeV * fm^-1) * fm^-1 * 1 = MeV  ✓
% T_{1/2} [s] = hbar [MeV·s] * ln(2) / Gamma_alpha [MeV]
Gamma  = 4 * (hbar_c^2 * ktilde^2 / (mu_c2 * k)) * uR_sq * chi^2;

if Gamma <= 0
    warning('tpa_halflife: Gamma <= 0（chi^2 = %.3e，uR_sq = %.3e）。', chi^2, uR_sq);
    T_half = Inf;
    return;
end

T_half = hbar_MeVs * log(2) / Gamma;
end
