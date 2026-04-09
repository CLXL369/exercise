function lam = find_lambda(r, VN0, VC, Q_c, G, L, hb2_2mu)
% FIND_LAMBDA  Find renormalisation factor lambda via Schrödinger eigenvalue condition.
%
%   lam = find_lambda(r, VN0, VC, Q_c, G, L, hb2_2mu)
%
% =========================================================================
% METHOD — Schrödinger eigenvalue condition (Bai & Ren 2018, Sec. 2)
% =========================================================================
%
%   The paper states: "the renormalisation factor lambda is determined by
%   requiring the alpha-core effective potential to reproduce the
%   experimentally measured Q_alpha value for L=0."
%
%   This means: find lambda such that the radial Schrödinger equation
%
%       [ -hbar^2/(2*mu) * d^2/dr^2  +  V_eff(r) ] u(r)  =  Q_c * u(r)
%
%   has a normalizable solution with exactly  N_target = (G-L)/2 = 8
%   interior nodes in the inner classically allowed well [0, r1].
%   This is the Wildermuth quantisation condition solved exactly by
%   numerical integration (not in the WKB approximation).
%
%   Algorithm — node-count bisection (Sturm-Liouville theorem):
%     1. Fix energy E = Q_c (the central experimental value).
%     2. Integrate the Schrödinger equation from r=0 outward by Numerov.
%     3. Count nodes of u(r) in [0, r1]  (the inner classically allowed region).
%     4. Increasing lambda deepens the well → more nodes fit at energy E=Q_c.
%     5. The sign of u just before r1 flips each time a new node enters the
%        well, enabling precise bisection to pin down the exact lambda*.
%
%   Inputs
%     r        – radial grid (fm), column vector starting at dr
%     VN0      – V_N(r)/lambda: unscaled nuclear potential (MeV)
%     VC       – V_C(r): Coulomb potential, independent of lambda (MeV)
%     Q_c      – central experimental Q_alpha (MeV)
%     G        – Wildermuth global quantum number (= 16)
%     L        – orbital angular momentum (= 0 for ground-state transitions)
%     hb2_2mu  – hbar^2/(2*mu) in MeV.fm^2
%   Output
%     lam – lambda such that E_{N_target}(lam) = Q_c  (Schrödinger condition)

dr       = r(2) - r(1);
N_target = (G - L) / 2;    % number of nodes = 8 for G=16, L=0   (eq. 18)

% ---- Coarse scan to locate the bracket [lam_lo, lam_hi] ----------------
% For each lambda, count nodes of u(r; E=Q_c) in the inner classically
% allowed well.  The transition from N_target-1 to N_target nodes marks
% the bracket that contains lambda*.
lam_arr  = linspace(0.25, 1.25, 60);
node_arr = zeros(1, numel(lam_arr));
u_arr    = zeros(1, numel(lam_arr));   % u at last well point (only sign used)

for i = 1:numel(lam_arr)
    V_eff = lam_arr(i) * VN0 + VC;
    if L > 0
        V_eff = V_eff + hb2_2mu * L*(L+1) ./ r.^2;
    end
    [node_arr(i), u_arr(i)] = nodes_and_u(r, V_eff, Q_c, hb2_2mu, dr);
end

% First index in scan where node count reaches N_target
idx_hi = find(node_arr >= N_target, 1);
if isempty(idx_hi)
    lam = lam_arr(end);
    warning('find_lambda: N_target=%d nodes not reached; returning lam=%.4f', N_target, lam);
    return;
end
if idx_hi == 1
    lam = lam_arr(1);
    warning('find_lambda: bracket not found below lam=%.3f; returning lam=%.4f', lam_arr(1), lam);
    return;
end

lam_lo = lam_arr(idx_hi - 1);
lam_hi = lam_arr(idx_hi);
u_lo   = u_arr(idx_hi - 1);

% ---- Bisection on sign of u(last well point) ----------------------------
% At lam_lo: N_target-1 nodes → sign is u_lo
% At lam_hi: N_target   nodes → sign is opposite to u_lo
% Bisect until |lam_hi - lam_lo| < 1e-8
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
% Count nodes of the Numerov solution u(r; E) in the classically allowed
% inner region [0, r1], and return u at the last grid point before r1.
%
% Schrödinger equation in reduced radial form:
%   u''(r) + g(r)*u(r) = 0,   g(r) = (E - V_eff(r)) / hb2_2mu
%   u(0) = 0,  u(dr) = dr   (regular at origin for L=0)
%
% Numerov recursion (4th-order, step h = dr):
%   u_{n+1} = [ 2*u_n*(1 - 5h^2/12 * g_n) - u_{n-1}*(1 + h^2/12 * g_{n-1}) ]
%             / (1 + h^2/12 * g_{n+1})

N = length(r);
g = (E - V_eff) / hb2_2mu;    % g(r): > 0 in well (classically allowed)

% Inner turning point: first index where g becomes non-positive
idx_r1 = find(g <= 0, 1);
if isempty(idx_r1) || idx_r1 <= 2
    n = 0;  u_last = r(1);  return;
end
n_steps = idx_r1 - 1;          % integrate r(1) through r(n_steps)

% g at the notional r=0 point: V_eff(0) ~ 0 for folded potentials
g0     = E / hb2_2mu;
u_prev = 0;                     % u(r = 0)
u_curr = r(1);                  % u(r = dr) = dr  [regular L=0 initial condition]
n      = 0;
s_prev = +1;                    % sign(u_curr) > 0

for j = 1 : n_steps - 1
    % g at previous, current, and next grid points
    g_prev = (j == 1) * g0 + (j > 1) * g(j - 1);
    g_curr = g(j);
    g_next = g(j + 1);

    denom  = 1.0 + (dr*dr / 12.0) * g_next;
    u_next = (2.0 * u_curr * (1.0 - (5.0/12.0)*dr*dr*g_curr) ...
              - u_prev * (1.0 + (dr*dr/12.0)*g_prev)) / denom;

    % Count sign changes (= nodes)
    s_next = sign(u_next);
    if s_next ~= 0 && s_next * s_prev < 0
        n      = n + 1;
        s_prev = s_next;
    end

    % Prevent floating-point overflow by renormalising every step
    % (sign is preserved; only magnitude is scaled)
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

u_last = u_curr;    % value (sign) of u just before the inner turning point
end
