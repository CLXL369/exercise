%% compute_table1.m
%
% Reproduce Table 1 of Bai & Ren, Eur. Phys. J. A (2018) 54:220
% "α clustering slightly above 100Sn in the light of the new experimental
%  data on the superallowed α decay"
%
% Method: Density-Dependent Cluster Model (DDCM) + Two-Potential Approach (TPA)
% Decay chain: 108Xe -> 104Te + alpha -> 100Sn + alpha
%
% Key equations used
%   V_alpha-core(r) = V_N(r) + V_C(r)                                 (eq. 2)
%   V_N via double-folding with M3Y interaction (lambda renorm.)       (eq. 3,5)
%   V_C via double-folding with Coulomb interaction                    (eq. 4,6)
%   T_{1/2} = hbar*ln2/Gamma,  Gamma=(4hbar^2 k~^2/mu/k)|phi chi|^2  (eq. 7)
%   lambda found by Schrodinger eigenvalue condition (eq. 18 exactly)
%   P_alpha = 1                                                        (sect. 2)
%
% Computation flow (per system, following the paper's logic):
%   Step 1: Build double-folding potentials VN0 and VC using Q_alpha_centre.
%           (VN0 = V_N/lambda; VC = Coulomb folded potential)
%   Step 2: Find lambda by requiring that the Schrodinger equation for the
%           alpha+core system gives an eigenvalue exactly equal to Q_alpha_centre
%           with N_target = (G-L)/2 = 8 nodes (Wildermuth condition, eq. 18).
%           This fixes the potential depth for the system.
%   Step 3: Use the *same* fixed V = lambda*VN0 + VC for all three Q_alpha values
%           [lower, centre, upper].  Only Q_alpha changes in the TPA formula (eq. 7);
%           the potential does not change.
%   Step 4: T_{1/2} = hbar*ln2 / Gamma_alpha  for each Q_alpha.
%
% Physical units: fm (lengths), MeV (energies), s (time)
%
% Usage:  octave compute_table1.m   OR   matlab compute_table1.m

clear; clc;

%% =========================================================================
%% Physical constants
%% =========================================================================
hbar_c    = 197.3269804;   % MeV.fm    (hbar * c)
e2        = 1.43997;       % MeV.fm    (e^2 = e^2/(4*pi*eps0) in natural units)
amu       = 931.494;       % MeV/c^2   (atomic mass unit)
hbar_MeVs = 6.58212e-22;   % MeV.s     (hbar)

%% =========================================================================
%% Decay systems
%% =========================================================================
% System 1: 108Xe -> 104Te + alpha
sys(1).label  = '108Xe -> 104Te';
sys(1).A_a    = 4;   sys(1).Z_a = 2;     % alpha particle
sys(1).A_c    = 104; sys(1).Z_c = 52;    % daughter core: 104Te
sys(1).L      = 0;                        % orbital angular momentum (g.s.)
sys(1).Q_c    = 4.6;                      % central Q_alpha (MeV), Bai&Ren eq.(1)
sys(1).Q_range = [4.4, 4.6, 4.8];        % [lower, centre, upper]  Q_alpha (MeV)

% System 2: 104Te -> 100Sn + alpha
sys(2).label  = '104Te -> 100Sn';
sys(2).A_a    = 4;   sys(2).Z_a = 2;
sys(2).A_c    = 100; sys(2).Z_c = 50;    % daughter core: 100Sn
sys(2).L      = 0;
sys(2).Q_c    = 5.1;                      % central Q_alpha (MeV)
sys(2).Q_range = [4.9, 5.1, 5.3];

%% =========================================================================
%% Radial grid
%% =========================================================================
dr = 0.04;     % fm  (step size)
Nr = 1200;     % grid points  (r_max = 48 fm)
r  = (1:Nr)' * dr;   % r(1) = 0.04 fm, r(Nr) = 48 fm

%% =========================================================================
%% Wildermuth global quantum number (eq. 18)
%% =========================================================================
% For 104Te and 108Xe: four valence nucleons in 0g_{7/2} (n_i=0, l_i=4)
% G = sum(2*n_i + l_i) = 4*4 = 16  =>  N_target = (G-L)/2 = 8 for L=0
G_wild = 16;

%% =========================================================================
%% Output header
%% =========================================================================
fprintf('\n');
fprintf('======================================================\n');
fprintf('  Bai & Ren (2018) Table 1 Reproduction\n');
fprintf('  DDCM + TPA,  P_alpha=1,  G=16\n');
fprintf('  lambda from Schrodinger eigenvalue condition\n');
fprintf('======================================================\n\n');
fprintf('%-22s  %8s  %8s  %14s\n', 'System', 'Q_a(MeV)', 'lambda', 'T_half(ns)');
fprintf('%s\n', repmat('-', 58, 1));

%% =========================================================================
%% Main loop over decay systems
%% =========================================================================
for is = 1:2
    s = sys(is);

    %% Reduced mass
    mu_c2   = s.A_a * s.A_c / (s.A_a + s.A_c) * amu;   % MeV/c^2
    hb2_2mu = hbar_c^2 / (2 * mu_c2);                    % MeV.fm^2

    %% Fermi density parameters for core nucleus (eq. 16)
    c_f  = 1.07 * s.A_c^(1/3);   % half-density radius (fm)
    a_f  = 0.54;                  % surface diffuseness (fm)
    rho0 = compute_rho0(s.A_c, c_f, a_f);

    fprintf('\n  [%s]  mu*c^2=%.2f MeV,  c=%.3f fm,  a=%.2f fm,  rho0=%.5f fm^-3\n', ...
            s.label, mu_c2, c_f, a_f, rho0);

    %% Step 1: Double-folding potentials (computed ONCE with central Q_alpha)
    % Using the central Q value for the M3Y exchange term J_EX = 276*(0.005*Q/A - 1)
    % ensures the potential is fixed for all three Q calculations below.
    [VN0, VC] = double_folding(r, s.A_a, s.Z_a, s.A_c, s.Z_c, ...
                               rho0, c_f, a_f, s.Q_c, e2);

    %% Step 2: Find lambda from Schrodinger eigenvalue condition (ONCE, central Q)
    % Requires E_{N_target}(lambda) = Q_c  exactly.
    % N_target = (G-L)/2 = 8 nodes  (Wildermuth condition, eq. 18)
    lam = find_lambda(r, VN0, VC, s.Q_c, G_wild, s.L, hb2_2mu);

    %% Full effective potential (fixed for this system)
    V = lam * VN0 + VC;
    if s.L > 0
        V = V + hb2_2mu * s.L*(s.L+1) ./ r.^2;
    end

    fprintf('  lambda = %.6f  (Schrodinger condition, Q_c = %.2f MeV)\n', lam, s.Q_c);

    %% Step 3: Compute T_{1/2} for each Q_alpha using the SAME fixed V
    % Only Q_a changes between iterations — the potential V is unchanged.
    % This reproduces the [lower, centre, upper] range of Table 1.
    for iq = 1:3
        Q_a = s.Q_range(iq);

        %% Step 4: Half-life via TPA (eq. 7)
        T_s  = tpa_halflife(r, V, Q_a, s.Z_a, s.Z_c, ...
                            mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, s.L);
        T_ns = T_s * 1e9;   % convert to nanoseconds

        fprintf('  %-22s  %8.2f  %8.6f  %14.4e\n', s.label, Q_a, lam, T_ns);
    end
end

fprintf('\n');
fprintf('------------------------------------------------------\n');
fprintf('  Paper Table 1 targets:\n');
fprintf('  108Xe->104Te:  T_half^th = (4.1-213)x10^3 ns  (centre 28e3 ns at Q=4.6)\n');
fprintf('  104Te->100Sn:  T_half^th = 7-166 ns           (centre 32 ns  at Q=5.1)\n');
fprintf('======================================================\n\n');
