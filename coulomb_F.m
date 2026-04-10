function F = coulomb_F(L, eta, rho)
% COULOMB_F  正则库仑函数 F_L(eta, rho)。
%
%   F = coulomb_F(L, eta, rho)
%
%   在经典禁戒区（rho < rho_tp）使用 WKB 近似，
%   在经典允许区（rho >= rho_tp）使用渐近正弦形式。
%   对 Sommerfeld 参数 eta >> 1 的情况（如本文研究的 alpha 衰变系统，eta ~ 10-15），
%   WKB 近似精度很高。
%
%   势垒内 WKB 公式（Blatt & Weisskopf；Mott & Massey）：
%       K(rho) = sqrt( 2*eta/rho - L*(L+1)/rho^2 - 1 )
%       I_C    = integral_{rho}^{rho_tp}  K(t)  dt
%       F_L(eta, rho) ≈ (1/2) * K(rho)^{-1/2} * exp(-I_C)
%
%   势垒外渐近形式：
%       theta_L = rho - eta*ln(2*rho) - L*pi/2 + sigma_L
%       F_L ~ sin(theta_L),  sigma_L = arg Gamma(L+1+i*eta)
%
%   输入参数
%     L    – 轨道角动量量子数（整数，>= 0）
%     eta  – Sommerfeld 参数（无量纲，> 0）
%     rho  – 无量纲变量 rho = k * R  （k 单位 fm^-1，R 单位 fm）
%   输出参数
%     F    – F_L(eta, rho) 的值

% 库仑经典转折点（近似，对 L=0 或小 L 精确）
% L(L+1)/rho^2 + 2*eta/rho - 1 = 0 => rho_tp ≈ eta + sqrt(eta^2 + L*(L+1))
rho_tp = eta + sqrt(eta^2 + L*(L+1));

if rho >= rho_tp
    %% 势垒外：渐近正弦形式
    sigma_L = imag(log(gamma(L + 1 + 1i*eta)));   % 库仑相移
    theta_L = rho - eta*log(2*rho) - L*pi/2 + sigma_L;
    F = sin(theta_L);
else
    %% 势垒内：WKB 近似
    K_rho = sqrt(max(2*eta/rho - L*(L+1)/(rho^2) - 1, 0));

    % 势垒积分 I_C = integral_rho^{rho_tp} K(t) dt
    Nt  = 4000;
    t   = linspace(rho, rho_tp, Nt);
    Kt  = sqrt(max(2*eta./t - L*(L+1)./t.^2 - 1, 0));
    I_C = trapz(t, Kt);

    % WKB：F ≈ (1/2) * K^{-1/2} * exp(-I_C)
    if K_rho > 0
        F = 0.5 * K_rho^(-0.5) * exp(-I_C);
    else
        % 在转折点处或非常靠近转折点；使用极限形式
        F = 0.5 * exp(-I_C);
    end
end
end
