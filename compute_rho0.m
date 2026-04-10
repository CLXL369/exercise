function rho0 = compute_rho0(A_c, c, a)
% COMPUTE_RHO0  计算 Fermi 核密度的归一化常数（公式14, 16）。
%
%   rho0 = compute_rho0(A_c, c, a)
%
%   求 rho0，使得：
%       4*pi * integral_0^inf  r^2 * rho0 / (1 + exp((r-c)/a))  dr  =  A_c
%
%   输入参数
%     A_c  – 子核质量数
%     c    – 半密度半径（fm）
%     a    – 表面弥散参数（fm）
%   输出参数
%     rho0 – 核心密度（fm^-3）

r   = linspace(0, c + 30*a, 10000)';
dr  = r(2) - r(1);
I   = 4*pi * sum( r.^2 ./ (1 + exp((r - c) ./ a)) ) * dr;
rho0 = A_c / I;
end
