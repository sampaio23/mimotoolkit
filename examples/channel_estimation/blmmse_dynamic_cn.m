tic;
channel_iterations = 2000;
noise_iterations = 1;
K = 2;
M = 4;
tau = K;
alpha = M;

SNRs = linspace(-30, 30, 21);

params.K = K;
params.M = M;
params.tau = tau;
params.channel_iterations = channel_iterations;
params.noise_iterations = noise_iterations;

params.alpha = alpha;
[~, Beff] = mtk_util_comparator_network('independent', params, 0);
[~, Beff_dynamic] = mtk_util_comparator_network('dynamic', params, 0);

params.channel_index = 1;
params.r_k = 0.5;
params.eta = rand(1,K);
H = mtk_generate_channel('kron-markov', params, 0);
Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params, 1);

params.C_h_real = 1/2*eye(2*M*K);
params.Phi = Phi;
params.quantizer = '1bit';

params_standard = params;
params_standard.Beff = Beff;
params_dynamic = params;
params_dynamic.Beff = Beff_dynamic;
clear params;

nmse_blmmse = zeros(size(SNRs));
nmse_blmmse_cn = zeros(size(SNRs));
nmse_blmmse_dynamic_cn = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_blmmse = 0;
    current_error_blmmse_cn = 0;
    current_error_blmmse_dynamic_cn = 0;

    params_standard.rho = mtk_util_db_to_linear(SNR);
    params_dynamic.rho = mtk_util_db_to_linear(SNR);
    params_blmmse = mtk_ce_prepare('blmmse', params_standard);
    params_blmmse_cn = mtk_ce_prepare('blmmse-cn', params_standard);
    params_blmmse_dynamic_cn = mtk_ce_prepare('blmmse-cn', params_dynamic);

    for channel_index=1:channel_iterations
        params_blmmse.H = H(:,:,channel_index);
        params_blmmse_cn.H = H(:,:,channel_index);
        params_blmmse_dynamic_cn.H = H(:,:,channel_index);
        params_blmmse.N = N(:,:,1);
        params_blmmse_cn.N = N(:,:,1);
        params_blmmse_dynamic_cn.N = N(:,:,1);

        H_hat = mtk_ce_estimate('blmmse', params_blmmse);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse = current_error_blmmse + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('blmmse-cn', params_blmmse_cn);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse_cn = current_error_blmmse_cn + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('blmmse-cn', params_blmmse_dynamic_cn);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse_dynamic_cn = current_error_blmmse_dynamic_cn + norm(E, 'fro')^2;
    end
    nmse_blmmse(snr_index) = current_error_blmmse/(channel_iterations*M*K);
    nmse_blmmse_cn(snr_index) = current_error_blmmse_cn/(channel_iterations*M*K);
    nmse_blmmse_dynamic_cn(snr_index) = current_error_blmmse_dynamic_cn/(channel_iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

plot(SNRs, nmse_blmmse, 'ro-', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, nmse_blmmse_cn, 'bo-', 'LineWidth', 1, 'MarkerSize', 8);
plot(SNRs, nmse_blmmse_dynamic_cn, 'ko-', 'LineWidth', 1, 'MarkerSize', 8);

plot(SNRs, (1-2/pi)*ones(length(SNRs)), '--')
plot(SNRs, (1-(5-2*sqrt(2))/pi)*ones(length(SNRs)), '--')

hold off;

legend('Numerical BLMMSE', 'Numerical BLMMSE-CN Independent', 'Numerical BLMMSE-CN Dynamic');
xlabel('\rho [dB]')
ylabel('MSE of channel estimation')
