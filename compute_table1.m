%% compute_table1.m
%
% 复现 Bai & Ren, Eur. Phys. J. A (2018) 54:220 中的表1
% "α clustering slightly above 100Sn in the light of the new experimental
%  data on the superallowed α decay"
%
% 方法：密度依赖团簇模型（DDCM）+ 双势方法（TPA）
% 衰变链：108Xe -> 104Te + alpha -> 100Sn + alpha
%
% 使用的关键方程
%   V_alpha-core(r) = V_N(r) + V_C(r)                                 （公式2）
%   V_N 通过 M3Y 相互作用的双折叠积分获得（lambda 归一化）             （公式3,5）
%   V_C 通过库仑相互作用的双折叠积分获得                               （公式4,6）
%   T_{1/2} = hbar*ln2/Gamma,  Gamma=(4hbar^2 k~^2/mu/k)|phi chi|^2  （公式7）
%   lambda 由薛定谔本征值条件精确求解（公式18）
%   P_alpha = 1                                                        （第2节）
%
% 计算流程（每个衰变系统，按论文逻辑）：
%   步骤1：用 Q_alpha 中心值计算双折叠势 VN0 和 VC。
%          （VN0 = V_N/lambda；VC = 库仑折叠势）
%   步骤2：通过要求 alpha+核薛定谔方程的本征值恰好等于 Q_alpha_centre，
%          且节点数 N_target = (G-L)/2 = 8（Wildermuth 条件，公式18），求 lambda。
%          此步骤固定系统的势阱深度。
%   步骤3：对三个 Q_alpha 值 [下限, 中心, 上限] 使用相同的固定势
%          V = lambda*VN0 + VC。TPA 公式（公式7）中只有 Q_alpha 变化，势不变。
%   步骤4：对每个 Q_alpha 计算 T_{1/2} = hbar*ln2 / Gamma_alpha。
%
% 物理单位：fm（长度），MeV（能量），s（时间）
%
% 使用方法：octave compute_table1.m  或  matlab compute_table1.m

clear; clc;

%% =========================================================================
%% 物理常数
%% =========================================================================
hbar_c    = 197.3269804;   % MeV·fm    （hbar * c）
e2        = 1.43997;       % MeV·fm    （e^2 = e^2/(4*pi*eps0)，自然单位）
amu       = 931.494;       % MeV/c^2   （原子质量单位）
hbar_MeVs = 6.58212e-22;   % MeV·s     （约化普朗克常数）

%% =========================================================================
%% 衰变系统参数
%% =========================================================================
% 系统1：108Xe -> 104Te + alpha
sys(1).label  = '108Xe -> 104Te';
sys(1).A_a    = 4;   sys(1).Z_a = 2;     % alpha 粒子
sys(1).A_c    = 104; sys(1).Z_c = 52;    % 子核：104Te
sys(1).L      = 0;                        % 轨道角动量（基态跃迁）
sys(1).Q_c    = 4.6;                      % Q_alpha 中心值（MeV），Bai&Ren 公式(1)
sys(1).Q_range = [4.4, 4.6, 4.8];        % [下限, 中心, 上限] Q_alpha（MeV）

% 系统2：104Te -> 100Sn + alpha
sys(2).label  = '104Te -> 100Sn';
sys(2).A_a    = 4;   sys(2).Z_a = 2;
sys(2).A_c    = 100; sys(2).Z_c = 50;    % 子核：100Sn
sys(2).L      = 0;
sys(2).Q_c    = 5.1;                      % Q_alpha 中心值（MeV）
sys(2).Q_range = [4.9, 5.1, 5.3];

%% =========================================================================
%% 径向网格
%% =========================================================================
dr = 0.04;     % fm  （步长）
Nr = 1200;     % 网格点数  （r_max = 48 fm）
r  = (1:Nr)' * dr;   % r(1) = 0.04 fm，r(Nr) = 48 fm

%% =========================================================================
%% Wildermuth 全局量子数（公式18）
%% =========================================================================
% 104Te 和 108Xe：四个价核子处于 0g_{7/2} 轨道（n_i=0, l_i=4）
% G = sum(2*n_i + l_i) = 4*4 = 16  =>  N_target = (G-L)/2 = 8（L=0）
G_wild = 16;

%% =========================================================================
%% 输出表头
%% =========================================================================
fprintf('\n');
fprintf('======================================================\n');
fprintf('  Bai & Ren (2018) 表1 复现\n');
fprintf('  DDCM + TPA,  P_alpha=1,  G=16\n');
fprintf('  lambda 由薛定谔本征值条件确定\n');
fprintf('======================================================\n\n');
fprintf('%-22s  %8s  %8s  %14s\n', '系统', 'Q_a(MeV)', 'lambda', 'T_half(ns)');
fprintf('%s\n', repmat('-', 58, 1));

%% =========================================================================
%% 主循环：遍历两个衰变系统
%% =========================================================================
for is = 1:2
    s = sys(is);

    %% 约化质量
    mu_c2   = s.A_a * s.A_c / (s.A_a + s.A_c) * amu;   % MeV/c^2
    hb2_2mu = hbar_c^2 / (2 * mu_c2);                    % MeV·fm^2

    %% 子核 Fermi 密度参数（公式16）
    c_f  = 1.07 * s.A_c^(1/3);   % 半密度半径（fm）
    a_f  = 0.54;                  % 表面弥散参数（fm）
    rho0 = compute_rho0(s.A_c, c_f, a_f);

    fprintf('\n  [%s]  mu*c^2=%.2f MeV,  c=%.3f fm,  a=%.2f fm,  rho0=%.5f fm^-3\n', ...
            s.label, mu_c2, c_f, a_f, rho0);

    %% 步骤1：双折叠势（用 Q_alpha 中心值计算一次）
    % 用中心 Q 值计算 M3Y 交换项 J_EX = 276*(0.005*Q/A - 1)，
    % 确保势在三个 Q 计算中保持不变。
    [VN0, VC] = double_folding(r, s.A_a, s.Z_a, s.A_c, s.Z_c, ...
                               rho0, c_f, a_f, s.Q_c, e2);

    %% 步骤2：由薛定谔本征值条件求 lambda（仅用中心 Q 一次）
    % 要求 E_{N_target}(lambda) = Q_c 精确成立。
    % N_target = (G-L)/2 = 8 个节点（Wildermuth 条件，公式18）
    lam = find_lambda(r, VN0, VC, s.Q_c, G_wild, s.L, hb2_2mu);

    %% 完整有效势（对该系统固定不变）
    V = lam * VN0 + VC;
    if s.L > 0
        V = V + hb2_2mu * s.L*(s.L+1) ./ r.^2;
    end

    fprintf('  lambda = %.6f  （薛定谔条件，Q_c = %.2f MeV）\n', lam, s.Q_c);

    %% 步骤3：用相同固定势 V 计算三个 Q_alpha 对应的半衰期
    % 三次迭代中只有 Q_a 变化，势 V 不变。
    % 这对应表1中的 [下限, 中心, 上限] 范围。
    for iq = 1:3
        Q_a = s.Q_range(iq);

        %% 步骤4：用 TPA 公式（公式7）计算半衰期
        T_s  = tpa_halflife(r, V, Q_a, s.Z_a, s.Z_c, ...
                            mu_c2, hb2_2mu, e2, hbar_c, hbar_MeVs, s.L);
        T_ns = T_s * 1e9;   % 转换为纳秒

        fprintf('  %-22s  %8.2f  %8.6f  %14.4e\n', s.label, Q_a, lam, T_ns);
    end
end

fprintf('\n');
fprintf('------------------------------------------------------\n');
fprintf('  论文表1 参考值：\n');
fprintf('  108Xe->104Te:  T_half^th = (4.1-213)x10^3 ns  (Q=4.6 时中心值 28e3 ns)\n');
fprintf('  104Te->100Sn:  T_half^th = 7-166 ns           (Q=5.1 时中心值 32 ns)\n');
fprintf('======================================================\n\n');
