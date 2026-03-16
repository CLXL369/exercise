function rho0 = compute_rho0(A_c, c, a)
% COMPUTE_RHO0  Normalisation constant for Fermi nuclear density (eq. 14, 16).
%
%   rho0 = compute_rho0(A_c, c, a)
%
%   Finds rho0 such that
%       4*pi * integral_0^inf  r^2 * rho0 / (1 + exp((r-c)/a))  dr  =  A_c
%
%   Inputs
%     A_c  – mass number of core nucleus
%     c    – half-density radius (fm)
%     a    – surface diffuseness (fm)
%   Output
%     rho0 – central density (fm^-3)

r   = linspace(0, c + 30*a, 10000)';
dr  = r(2) - r(1);
I   = 4*pi * sum( r.^2 ./ (1 + exp((r - c) ./ a)) ) * dr;
rho0 = A_c / I;
end
