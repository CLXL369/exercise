function lam = find_lambda(r, VN0, VC, Q_a, G, L, hb2_2mu)
% FIND_LAMBDA  Find renormalisation factor lambda via WKB Wildermuth condition.
%
%   lam = find_lambda(r, VN0, VC, Q_a, G, L, hb2_2mu)
%
%   Requires the WKB quantisation condition for the inner classically-allowed
%   well (eq. 18 combined with Bohr-Sommerfeld with hard-wall at r=0):
%
%       integral_0^{r1}  k(r)  dr  =  (G/2 + 3/4) * pi
%
%   where  k(r) = sqrt( (Q_a - V_eff(r)) / hb2_2mu )
%   and    r1   is the first inner classical turning point (V_eff = Q_a).
%   G = 16 gives target = 8.75*pi  (N=8 nodes, L=0).
%
%   Lambda is found by bisection on [lam_lo, lam_hi].
%
%   Inputs
%     r        – radial grid (fm), column vector
%     VN0      – V_N(r)/lambda, nuclear potential without lambda (MeV)
%     VC       – V_C(r), Coulomb potential (MeV)
%     Q_a      – decay Q-value (MeV)
%     G        – Wildermuth global quantum number (= 16)
%     L        – orbital angular momentum (= 0 for ground-state transitions)
%     hb2_2mu  – hbar^2/(2*mu) (MeV.fm^2)
%   Output
%     lam – lambda such that WKB phase = (G/2+3/4)*pi

dr     = r(2) - r(1);
target = (G/2 + 3/4) * pi;   % = 8.75*pi for G=16, L=0

% Bisection bounds: lambda in [0.5, 1.2] (paper states lambda ~ 0.5-0.9)
lam_lo = 0.50;
lam_hi = 1.20;

% Verify bracketing
I_lo = wkb_phase_inner(r, lam_lo*VN0 + VC, Q_a, L, hb2_2mu, dr);
I_hi = wkb_phase_inner(r, lam_hi*VN0 + VC, Q_a, L, hb2_2mu, dr);

% If bracket fails, extend search range
if I_lo >= target
    lam_lo = 0.20;
end
if I_hi <= target
    lam_hi = 1.60;
end

% Bisection (100 iterations -> precision ~1e-7 in lambda)
for iter = 1:100
    lam = 0.5 * (lam_lo + lam_hi);
    I   = wkb_phase_inner(r, lam*VN0 + VC, Q_a, L, hb2_2mu, dr);
    if I < target
        lam_lo = lam;    % phase too small -> need deeper well -> larger lambda
    else
        lam_hi = lam;
    end
    if (lam_hi - lam_lo) < 1e-8
        break;
    end
end
lam = 0.5 * (lam_lo + lam_hi);
end

% -------------------------------------------------------------------------
function I = wkb_phase_inner(r, V, Q_a, L, hb2_2mu, dr)
% Compute WKB phase integral from 0 to the first inner turning point r1.
%   I = integral_0^{r1}  sqrt( max(Q_a - V_eff(r), 0) / hb2_2mu )  dr
V_eff = V;
if L > 0
    V_eff = V + hb2_2mu * L*(L+1) ./ r.^2;
end

k2 = (Q_a - V_eff) / hb2_2mu;    % k^2(r); positive in classically allowed

% First index where the well becomes classically forbidden (k^2 <= 0)
idx_turn = find(k2 <= 0, 1);
if isempty(idx_turn)
    idx_turn = length(r) + 1;     % never forbidden -> integrate to end
end

% Integrate from r(1) to r(idx_turn - 1)
if idx_turn <= 2
    I = 0;
    return;
end
n_int = idx_turn - 1;
k_int = sqrt(max(k2(1:n_int), 0));
I     = sum(k_int) * dr;
end
