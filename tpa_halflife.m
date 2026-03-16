function T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
% TPA_HALFLIFE  Compute alpha-decay half-life via Two-Potential Approach (eq. 7).
%
%   T_half = tpa_halflife(r, V, Q_a, Z_a, Z_c, mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, L)
%
%   TPA formula (Bai & Ren eq. 7, P_alpha = 1):
%       Gamma_alpha = (4*hbar^2*ktilde^2) / (mu*k)  *  |phi_L(R)|^2  *  |chi_L(kR)|^2
%       T_{1/2}     = hbar * ln2 / Gamma_alpha
%
%   Inner wave function at R computed using the WKB approximation (valid for
%   eta >> 1, i.e., long-lived states where the barrier is wide):
%
%       |phi_L(R)|^2 = exp(-2*G_inner) / (2 * ktilde * I_inv)
%
%   where
%     G_inner = integral_{r1}^{R}  kappa(r)  dr      inner barrier action
%     I_inv   = integral_0^{r1}    k(r)^{-1} dr      WKB normalisation integral
%     k(r)    = sqrt((Q_a - V(r)) / hb2_2mu)         inner well wave number
%     kappa(r)= sqrt((V(r) - Q_a) / hb2_2mu)         inner barrier decay rate
%     r1      = first classical turning point (V(r1) = Q_a)
%     R       = separation radius  (barrier maximum, V'(R) = 0)
%     k       = sqrt(2*mu*Q_a) / hbar                asymptotic wave number
%     ktilde  = sqrt(2*mu*(V(R)-Q_a)) / hbar         imaginary wave number at R
%     chi_L   = F_L(eta, k*R)   regular Coulomb function (WKB, eq. 6)
%     eta     = Z_a*Z_c*e^2*mu / (hbar^2*k)          Sommerfeld parameter
%
%   The direct outward Numerov integration picks up the exponentially growing
%   WKB mode in the inner barrier [r1, R] and gives an overestimate of u(R)
%   by a factor ~ exp(G_inner). The WKB formula corrects this by construction.
%
%   All lengths in fm, energies in MeV.

dr = r(2) - r(1);

%% ---- Locate separation radius R (barrier maximum, V'(R)=0) -----------
dV      = diff(V) / dr;
r_m     = 0.5*(r(1:end-1) + r(2:end));
mask    = r_m > 5.0;
dV_m    = dV(mask);
idx_mask = find(mask);

sgn_chg = find(dV_m(1:end-1) > 0 & dV_m(2:end) < 0, 1, 'last');
if ~isempty(sgn_chg)
    idx_R = idx_mask(sgn_chg);
else
    sub_idx = find(r > 5.0);
    [~, loc] = max(V(sub_idx));
    idx_R = sub_idx(loc);
end
R = r(idx_R);

if V(idx_R) <= Q_a
    warning('tpa_halflife: V(R) <= Q_a; barrier not found.');
    T_half = Inf;
    return;
end

%% ---- Wave numbers at R -----------------------------------------------
k      = sqrt(2 * mu_c2 * Q_a)              / hbar_c;   % fm^-1
ktilde = sqrt(2 * mu_c2 * (V(idx_R) - Q_a)) / hbar_c;   % fm^-1

%% ---- Inner classical turning point r1 (first crossing V = Q_a) ------
idx_r1 = find(V(1:idx_R) > Q_a, 1);   % first grid point where V > Q_a
if isempty(idx_r1) || idx_r1 <= 1
    warning('tpa_halflife: inner turning point not found.');
    T_half = Inf;
    return;
end

%% ---- WKB normalisation integral  I_inv = int_0^{r1} k(r)^{-1} dr ---
% k(r) = sqrt((Q_a - V(r)) / hb2_2mu)  in the classically allowed inner well
k2_in  = (Q_a - V(1:idx_r1-1)) / hb2_2mu;   % k^2 > 0 for r < r1
k_in   = sqrt(max(k2_in, 0));
k_safe = max(k_in, 1e-6);                     % floor avoids 1/0 near r1

% Analytic endpoint correction for integrable singularity at r1:
% k(r) ~ sqrt(V'(r1) * (r1-r) / hb2_2mu) near r1
% int_{r1-dr}^{r1} k^{-1} dr ~ 2*sqrt(hb2_2mu / V'(r1)) * sqrt(dr)
V_prime  = (V(idx_r1) - V(idx_r1-1)) / dr;   % dV/dr at turning point (> 0)
if V_prime > 0
    endpoint = 2 * sqrt(hb2_2mu / V_prime) * sqrt(dr) ...
             - sqrt(hb2_2mu / (V_prime * dr));   % subtract double-counted last cell
else
    endpoint = 0;
end
I_inv = sum(1.0 ./ k_safe) * dr + endpoint;

%% ---- Inner barrier action  G_inner = int_{r1}^{R} kappa(r) dr ------
kap_sq  = (V(idx_r1:idx_R) - Q_a) / hb2_2mu;
kap     = sqrt(max(kap_sq, 0));
G_inner = sum(kap) * dr;

%% ---- WKB inner wave function amplitude at R -------------------------
% |phi_L(R)|^2 = exp(-2*G_inner) / (2 * ktilde * I_inv)
uR_sq = exp(-2 * G_inner) / (2 * ktilde * I_inv);

%% ---- Sommerfeld parameter and Coulomb function -----------------------
eta = Z_a * Z_c * e2 * mu_c2 / (hbar_c^2 * k);
chi = coulomb_F(L, eta, k * R);          % F_L(eta, k*R)

%% ---- Decay width and half-life (eq. 7) --------------------------------
% Gamma [MeV] = (4*hbar_c^2*ktilde^2/(mu_c2*k)) * |phi(R)|^2 * chi^2
% Units: MeV^2.fm^2.fm^-2 / (MeV.fm^-1) * fm^-1 * 1 = MeV  ✓
Gamma  = 4 * (hbar_c^2 * ktilde^2 / (mu_c2 * k)) * uR_sq * chi^2;
T_half = hbar_MeVs * log(2) / Gamma;
end
