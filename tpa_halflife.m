function T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
% TPA_HALFLIFE  Alpha-decay half-life via the Two-Potential Approach (TPA).
%
%   T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
%
% =========================================================================
% THEORY — Two-Potential Approach (TPA), Bai & Ren (2018) Sec. 2
% =========================================================================
%
% The total alpha–core potential is split into two regions at the channel
% radius R (the outer maximum of V):
%
%   Region I  (inner, r <= R): V1(r) = V_N(r) + V_C(r)   [full potential]
%   Region II (outer, r > R ): V2(r) = V_C(r)             [Coulomb only]
%
% ---- Eq. (7): decay width -----------------------------------------------
%   Gamma_alpha = (4 * hbar^2 * ktilde^2) / (mu * k)
%                 * |phi_L(R)|^2
%                 * |chi_L(eta, kR)|^2
%
%   where (all in fm, MeV units):
%     k      = sqrt(2*mu*Q_alpha) / hbar          asymptotic wave number
%     ktilde = sqrt(2*mu*(V(R)-Q_alpha)) / hbar   imaginary wave number at R
%     phi_L(R) = inner wave function at R, normalised integral_0^R |phi|^2 dr = 1
%     chi_L  = F_L(eta, kR) = regular Coulomb wave function [eq. (6) region II]
%     eta    = Z_alpha * Z_core * e^2 * mu / (hbar^2 * k)   [Sommerfeld param]
%
%   Half-life: T_{1/2} = hbar * ln(2) / Gamma_alpha        [P_alpha = 1]
%
% ---- |phi_L(R)|^2: WKB inner wave function at R -------------------------
%   A direct outward Numerov integration from r=0 is NOT used for |phi(R)|
%   because it captures the exponentially *growing* WKB mode in the inner
%   barrier [r1, R], over-estimating |phi(R)| by exp(G_inner).
%
%   Instead, the standard TPA implementation uses the WKB connection formula
%   to match the inner-well solution through the inner barrier:
%
%       |phi_L(R)|^2 = exp(-2 * G_inner) / (2 * ktilde * I_inv)
%
%   where:
%     r1      = inner classical turning point: first r where V(r) = Q_alpha
%     G_inner = integral_{r1}^{R}  kappa(r)  dr    [inner barrier action]
%     kappa(r)= sqrt((V(r)-Q_alpha)/hb2_2mu)        [local decay rate in barrier]
%     I_inv   = integral_0^{r1}  1/k(r)  dr         [WKB norm. integral]
%     k(r)    = sqrt((Q_alpha-V(r))/hb2_2mu)         [inner well wave number]
%
%   This is equivalent to eq. (7) with phi_L computed via WKB, which is
%   accurate because the Sommerfeld parameter eta >> 1 in these systems.
%
% =========================================================================
%   All lengths in fm, energies in MeV.

dr = r(2) - r(1);

%% ---- Step 1: Locate separation radius R --------------------------------
% R is the outermost maximum of V(r) for r > 5 fm (eq. 7 channel radius).
% We detect the last sign change + -> - in dV/dr beyond r=5 fm.
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
    warning('tpa_halflife: V(R) <= Q_a; barrier top not found.');
    T_half = Inf;
    return;
end

%% ---- Step 2: Wave numbers at R (eq. 7) ---------------------------------
% k      = sqrt(2*mu*Q_a)          / hbar    [asymptotic, region II]
% ktilde = sqrt(2*mu*(V(R)-Q_a))   / hbar    [imaginary, at barrier top]
k      = sqrt(2 * mu_c2 * Q_a)              / hbar_c;
ktilde = sqrt(2 * mu_c2 * (V(idx_R) - Q_a)) / hbar_c;

%% ---- Step 3: Inner classical turning point r1 --------------------------
% r1 is the first radial point where V(r) rises above Q_alpha (going outward),
% i.e., the first root of V(r) = Q_alpha in [0, R].
idx_r1 = find(V(1:idx_R) > Q_a, 1);
if isempty(idx_r1) || idx_r1 <= 1
    warning('tpa_halflife: inner turning point not found.');
    T_half = Inf;
    return;
end

%% ---- Step 4: WKB normalisation integral I_inv = int_0^{r1} k^{-1} dr --
% k(r) = sqrt((Q_a - V(r)) / hb2_2mu)  in the classically allowed well.
% The integrand 1/k(r) is bounded (k >= ~0.25 fm^-1 in these systems), so
% a plain rectangular sum gives a reliable result without singularity treatment.
k2_in  = (Q_a - V(1:idx_r1-1)) / hb2_2mu;   % positive inside the well
k_in   = sqrt(max(k2_in, 0));
k_safe = max(k_in, 1e-6);                     % safety floor (never triggered for typical Q)
I_inv  = sum(1.0 ./ k_safe) * dr;             % fm^2 (units: [fm^-1]^-1 * fm = fm^2 ✓)

if I_inv <= 0
    warning('tpa_halflife: I_inv <= 0; unable to normalise inner wave function.');
    T_half = Inf;
    return;
end

%% ---- Step 5: Inner barrier action G_inner = int_{r1}^{R} kappa dr -----
% kappa(r) = sqrt((V(r)-Q_a)/hb2_2mu)  in the classically forbidden region.
kap_sq  = (V(idx_r1:idx_R) - Q_a) / hb2_2mu;
kap     = sqrt(max(kap_sq, 0));
G_inner = sum(kap) * dr;

%% ---- Step 6: WKB inner wave function amplitude at R --------------------
% |phi_L(R)|^2 = exp(-2*G_inner) / (2 * ktilde * I_inv)   [always positive]
uR_sq = exp(-2 * G_inner) / (2 * ktilde * I_inv);

%% ---- Step 7: Sommerfeld parameter and Coulomb function (eq. 6) --------
% eta = Z_a*Z_c*e^2*mu / (hbar^2*k)  [dimensionless Sommerfeld parameter]
% chi = F_L(eta, kR)                 [regular Coulomb function, region II]
eta = Z_a * Z_c * e2 * mu_c2 / (hbar_c^2 * k);
chi = coulomb_F(L, eta, k * R);

%% ---- Step 8: Decay width and half-life (eq. 7) -------------------------
% Gamma_alpha [MeV] = 4*(hbar^2*ktilde^2/(mu*k)) * |phi(R)|^2 * |chi|^2
%   Unit check: MeV^2.fm^2 * fm^-2 / (MeV * fm^-1) * fm^-1 * 1 = MeV  ✓
% T_{1/2} [s] = hbar [MeV.s] * ln(2) / Gamma_alpha [MeV]
Gamma  = 4 * (hbar_c^2 * ktilde^2 / (mu_c2 * k)) * uR_sq * chi^2;

if Gamma <= 0
    warning('tpa_halflife: Gamma <= 0 (chi^2 = %.3e, uR_sq = %.3e).', chi^2, uR_sq);
    T_half = Inf;
    return;
end

T_half = hbar_MeVs * log(2) / Gamma;
end
