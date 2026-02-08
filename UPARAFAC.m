% 基于酉 PARAFAC 的联合 DOA 与载频估计算法（ULA 场景）
%
% 本代码为论文方法的复现实现，用于学习与科研验证。
%
% 【参考文献】
% Lingyun Xu, Fangqing Wen, Xiaofei Zhang,
% “A Novel Unitary PARAFAC Algorithm for Joint DOA and
% Frequency Estimation”,
% IEEE Communications Letters,
% vol. 23, no. 4, pp. xxx–xxx, Apr. 2019.
% DOI: 10.1109/LCOMM.2019.2896593
%
% 【声明】
% 本代码为作者基于公开论文的独立复现版本，
% 并非论文作者提供的官方代码。
% 在仿真参数设置、工具函数实现及部分细节上
% 可能与原文或作者实现存在差异。
%
% 【方法说明】
% - 采用前向–后向扩展（Forward-Backward Averaging）
% - 结合酉变换（Unitary Transform），将复值模型
%   转换为实值 PARAFAC 张量模型
% - 通过三线性分解实现 DOA 与载频的联合估计
%
% 【可调参数】
% M, N, tau, d, f, theta, L, SNR 等参数
% 可根据实验需求自行修改
%
% 【用途】
% 仅用于学习、研究与算法复现实验




clc; clear all; close all;
%% 信号模型参数
M = 12;                     % 阵元数
N = 10;                    % 延迟通道数
d = 15;                   % 阵元间距，不是越小越好，在半波长左右最好
c = 3e8;                   % 光速
tau = 0.01e-6;             % 延迟时间
theta = [-5, 5, 15];      % DOA角度（度）
f = [5.5e6, 8.5e6, 9.5e6]; % 信源频率（Hz）；在满足条件的情况下，f越高越准，包括最小的f
K = length(theta);
L = 100;                   % 快拍数
SNR = 10;                  % 信噪比（dB）
item = 50;

%% ========================如何调参数=============================
%波长≥2d，满足不等式条件，可调参数τ、N、f
%波长约束：λ = c/f ≥ 1米 ⇒ f ≤ 3e8 Hz；不等式约束：(N-1)·τ < 1/max(f)
%参数范围：f ≤ 300MHz，τ很小(1e-9量级)，N适中(6-12)
%先确定τ：选择较小的τ值(1e-9 ~ 5e-10)
%再确定N：根据τ选择适当的N(6-10)
%最后调f：确保f ≤ 1/[(N-1)·τ] 且 f ≤ 300MHz
% 计算波长
lambda = c ./ f;          

% 验证不等式条件 0 < (N-1)*tau < 1/max(fk)
left_bound = 0;
mid_value = (N-1) * tau;
right_bound = 1 / max(f);
fprintf('=== 不等式条件验证 ===\n');
fprintf('条件: 0 < (N-1)*tau < 1/max(fk)\n');
fprintf('计算: 0 < (%d-1)*%.2e < 1/%.1e\n', N, tau, max(f));
fprintf('结果: 0 < %.2e < %.2e\n', mid_value, right_bound);
if left_bound < mid_value && mid_value < right_bound
    fprintf('✓ 条件满足！\n');
else
    fprintf('✗ 条件不满足！\n');
end
fprintf('\n参数汇总:\n');
fprintf('阵元数 M = %d\n', M);
fprintf('延迟通道数 N = %d\n', N);
fprintf('延迟时间 tau = %.2e s\n', tau);
fprintf('频率范围: %.1f - %.1f MHz\n', min(f)/1e6, max(f)/1e6);
fprintf('波长范围: %.1f - %.1f 米\n', min(lambda), max(lambda));
%%===============================================================



%% 构建X_bar
S = randn(K,L) + 1j*randn(K,L); 
alpha = 2*pi*d*f.*sind(theta)/c;
beta = 2*pi*f*tau; 
As = exp(-1j * (0:M-1)' * alpha);                         % 非对称空间导向矩阵 M×K
Af = exp(-1j * (0:N-1)' * beta);                          % 非对称频率导向矩阵 N×K
Asf = khatriRao(Af, As); 
Gamma_MN = fliplr(eye(M*N));
Gamma_L = fliplr(eye(L));
for item_num = 1:item
	X = awgn(Asf*S,SNR,'measured','dB');
	X_bar = [X,Gamma_MN*conj(X)*Gamma_L];

	%% 实值张量构建与PARAFAC
	U_M = build_unitary_matrix(M);                          % M×M
	U_N = build_unitary_matrix(N);                          % N×N
	U_2L = build_unitary_matrix(2*L);                       % 2L×2L
	Y_real = kron(U_N', U_M') * X_bar * U_2L;  % 理论上是实数
	%Y_real = real(Y_real);                    % 清除浮点虚部误差（非常重要）
	Y_tensor = reshape(Y_real, M, N, 2*L);    % 转为三阶张量
	[E_Asr, E_Afr, E_Sr, ~] = comfac(Y_tensor, K);  % 三线性分解


	%% 使用 phase() 提取相位差（逐列计算）
	% 酉逆变换（还原频率与空间方向）
	U_fk = U_N * E_Afr;                                     % 频率方向 N×K
	U_sk = U_M * E_Asr;                                     % 空间方向 M×K
	h_fk = zeros(N, K);
	h_sk = zeros(M, K);
	for k = 1:K
		h_fk(:,k) = -phase(U_fk(:,k) / U_fk(1,k));
		h_sk(:,k) = -phase(U_sk(:,k) / U_sk(1,k));
	end


	%% 构造频率估计设计矩阵，最小二乘估计
	u = (0:N-1)' * (2*pi*tau);
	P1 = [ones(N,1),u];
	c1 = pinv(P1) * h_fk;                                    % N×K → 2×K
	estimated_freq = c1(2,:);                  % 提取频率估计值

	%% 构造角度估计设计矩阵，最小二乘估计
	for k = 1:K
		v = (0:M-1)' * (estimated_freq(k)/c);     % 构造 v，无 sinθ
		P2 = [ones(M,1), v];             % 构造 P2 = [1, 2πd·v]
		c2 = pinv(P2) * h_sk(:,k);                % 最小二乘拟合
		estimated_theta(k) = asind(c2(2) / (2*pi*d));  % 解出 sinθ → θ
    end
    
    
	%%画图
    plot(estimated_freq,estimated_theta,'k*');hold on;
end

%%可视化
H(1)=plot(estimated_freq,estimated_theta,'k*');hold on; %只是为了拿句柄，与上一个plot重复了
H(2)=plot(f,theta,'rx','markersize',28);grid on;
xlabel('Frequency'),ylabel('DOA');
legend([H(1),H(2)],'Esimated','Ture')


%%======辅助函数：构建酉矩阵======
function U = build_unitary_matrix(k)
    if mod(k,2) == 0
        % 偶数维度
        k_half = k/2;
        I = eye(k_half);
        J = fliplr(eye(k_half));
        U = (1/sqrt(2)) * [I, 1j*I; J, -1j*J];
    else
        % 奇数维度  
        k_half = floor(k/2);
        I = eye(k_half);
        J = fliplr(eye(k_half));
        zero_col = zeros(k_half,1);
        U = (1/sqrt(2)) * [I, zero_col, 1j*I; 
                           zero_col', sqrt(2), zero_col';
                           J, zero_col, -1j*J];
    end
end