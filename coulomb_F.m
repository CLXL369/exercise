function F = coulomb_F(L, eta, rho)
% COULOMB_F  Regular Coulomb function F_L(eta, rho).
%
%   F = coulomb_F(L, eta, rho)
%
%   Uses the WKB approximation in the classically forbidden region (rho < rho_tp)
%   and the asymptotic sinusoidal form in the allowed region (rho >= rho_tp).
%   The WKB approximation is highly accurate for the Sommerfeld parameter eta >> 1,
%   which holds for the alpha-decay systems studied here (eta ~ 10-15).
%
%   WKB formula inside barrier (Blatt & Weisskopf; Mott & Massey):
%       K(rho) = sqrt( 2*eta/rho - L*(L+1)/rho^2 - 1 )
%       I_C    = integral_{rho}^{rho_tp}  K(t)  dt
%       F_L(eta, rho) ≈ (1/2) * K(rho)^{-1/2} * exp(-I_C)
%
%   Asymptotic form outside barrier:
%       theta_L = rho - eta*ln(2*rho) - L*pi/2 + sigma_L
%       F_L ~ sin(theta_L),  sigma_L = arg Gamma(L+1+i*eta)
%
%   Inputs
%     L    – orbital angular momentum quantum number (integer >= 0)
%     eta  – Sommerfeld parameter (dimensionless, > 0)
%     rho  – dimensionless argument  rho = k * R  (k in fm^-1, R in fm)
%   Output
%     F    – value of F_L(eta, rho)

% Classical turning point for Coulomb (approximate, valid for L=0 or small L)
% L(L+1)/rho^2 + 2*eta/rho - 1 = 0 => rho_tp ≈ eta + sqrt(eta^2 + L*(L+1))
rho_tp = eta + sqrt(eta^2 + L*(L+1));

if rho >= rho_tp
    %% Outside barrier: asymptotic sinusoidal form
    sigma_L = imag(log(gamma(L + 1 + 1i*eta)));   % Coulomb phase shift
    theta_L = rho - eta*log(2*rho) - L*pi/2 + sigma_L;
    F = sin(theta_L);
else
    %% Inside barrier: WKB approximation
    K_rho = sqrt(max(2*eta/rho - L*(L+1)/(rho^2) - 1, 0));

    % Barrier integral I_C = integral_rho^{rho_tp} K(t) dt
    Nt  = 4000;
    t   = linspace(rho, rho_tp, Nt);
    Kt  = sqrt(max(2*eta./t - L*(L+1)./t.^2 - 1, 0));
    I_C = trapz(t, Kt);

    % WKB: F ≈ (1/2) * K^{-1/2} * exp(-I_C)
    if K_rho > 0
        F = 0.5 * K_rho^(-0.5) * exp(-I_C);
    else
        % At or very near the turning point; use limiting form
        F = 0.5 * exp(-I_C);
    end
end
end
