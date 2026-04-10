function [VN0, VC] = double_folding(r, A_a, Z_a, A_c, Z_c, rho0, c_f, a_f, Q_a, e2)
% DOUBLE_FOLDING  计算核力（VN0）和库仑（VC）双折叠势。
%
%   [VN0, VC] = double_folding(r, A_a, Z_a, A_c, Z_c, rho0, c_f, a_f, Q_a, e2)
%
%   核力势（不含 lambda）：
%       V_N(r) = lambda * VN0(r)
%       V_N(r) = lambda * int dr_a dr_c rho_a(r_a) rho_c(r_c) v_n(Q_a, s)
%   库仑势：
%       V_C(r) = int dr_a dr_c rho~_a(r_a) rho~_c(r_c) v_c(s)
%
%   两者均通过傅里叶卷积定理计算（公式3-6）：
%       V(r) = 1/(2*pi^2*r) * int_0^inf  k * f_a(k) * f_c(k) * v_tilde(k) * sin(k*r) dk
%
%   密度分布
%     alpha 粒子：rho_a(r) = 0.422875 * exp(-0.7024*r^2)     [公式13，归一化到 A_a=4]
%     子核：      rho_c(r) = rho0 / (1 + exp((r-c)/a))       [公式14]
%
%   M3Y 核子-核子相互作用（公式5）：
%     v_n(s) = 7999*exp(-4s)/(4s) - 2134*exp(-2.5s)/(2.5s) + J_EX*delta(s)
%     J_EX   = 276*(0.005*Q_a/A_a - 1)   [MeV·fm^3]
%
%   库仑相互作用（公式6）：
%     v_c(s) = e^2/s，电荷密度 rho~ = (Z/A)*rho_nucleon
%
%   所有长度单位：fm；能量单位：MeV。

%% ---- k 网格 -----------------------------------------------------------
Nk    = 2000;
k_max = 12.0;            % fm^-1  （Yukawa + 高斯衰减使被积函数在此后可忽略）
k     = linspace(0.003, k_max, Nk)';
dk    = k(2) - k(1);

%% ---- Alpha 粒子密度形状因子（解析，高斯型）-------------------
% rho_a(r) = 0.422875 * exp(-0.7024*r^2)  =>  f_a(k) = A_a * exp(-k^2/(4*beta))
beta  = 0.7024;          % fm^-2
f_a   = A_a * exp(-k.^2 / (4*beta));

%% ---- 子核密度形状因子（Fermi 型，数值正弦变换）--------
% f_c(k) = (4*pi/k) * int_0^inf  r * rho_c(r) * sin(k*r)  dr
Nr_d  = 2000;
r_d   = linspace(0.005, c_f + 25*a_f, Nr_d)';
dr_d  = r_d(2) - r_d(1);
rho_c = rho0 ./ (1 + exp((r_d - c_f) ./ a_f));

% 矩阵乘积：sin_mat(i,j) = sin(k(i)*r_d(j))，大小 Nk x Nr_d
sin_mat = sin(k * r_d');                                    % Nk x Nr_d
f_c     = (4*pi ./ k) .* (sin_mat * (r_d .* rho_c) * dr_d); % Nk x 1

%% ---- M3Y 相互作用的傅里叶变换（公式5）----------------------
% Yukawa 项  v(s) = C*exp(-alpha*s)/s  =>  v_tilde(k) = 4*pi*C / (k^2 + alpha^2)
% v_n1: C1 = 7999/4  MeV·fm,  alpha1 = 4  fm^-1
% v_n2: C2 =-2134/2.5 MeV·fm, alpha2 = 2.5 fm^-1
vt_1  = 4*pi * ( 7999/4   ) ./ (k.^2 + 16.00);  % MeV·fm^3
vt_2  = 4*pi * (-2134/2.5 ) ./ (k.^2 +  6.25);  % MeV·fm^3
J_EX  = 276 * (0.005*Q_a/A_a - 1);               % MeV·fm^3（零程交换项）
vt_n  = vt_1 + vt_2 + J_EX;                       % 总 M3Y 傅里叶变换

%% ---- 库仑相互作用的傅里叶变换（公式6）-------------------
% v_c(s) = e^2/s  =>  v_tilde_c(k) = 4*pi*e^2 / k^2
% 电荷密度：rho~ = (Z/A)*rho_nucleon
vt_c   = 4*pi * e2 ./ k.^2;
f_a_ch = (Z_a / A_a) * f_a;    % alpha 电荷密度形状因子
f_c_ch = (Z_c / A_c) * f_c;    % 子核电荷密度形状因子

%% ---- 逆变换到实空间 ---------------------------------
% V(r) = dk/(2*pi^2) * sum_k  [k * F(k) * sin(k*r)] / r
Nr_r  = length(r);
r_row = r(:)';                          % 1 x Nr_r

% sin 矩阵：sin_kr(i,j) = sin(k(i)*r(j))，大小 Nk x Nr_r
sin_kr = sin(k * r_row);                 % Nk x Nr_r

FN  = k .* f_a    .* f_c    .* vt_n;    % 核力被积函数  (Nk x 1)
FC  = k .* f_a_ch .* f_c_ch .* vt_c;   % 库仑被积函数  (Nk x 1)

VN0 = (dk / (2*pi^2)) * (FN' * sin_kr) ./ r_row;   % 1 x Nr_r
VC  = (dk / (2*pi^2)) * (FC' * sin_kr) ./ r_row;   % 1 x Nr_r

VN0 = VN0(:);   % 列向量，单位 MeV
VC  = VC(:);    % 列向量，单位 MeV
end
