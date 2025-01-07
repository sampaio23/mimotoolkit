tic;
iterations = 4000;
K = 8;
M = 8;
tau = 8;
alpha = 8;

SNRs = linspace(-10, 20, 9);

params.K = K;
params.M = M;
params.tau = tau;
params.iterations = iterations;

params.alpha = alpha;
[params.B, params.Beff] = mtk_util_comparator_network('independent', params);

H = mtk_generate_channel('rayleigh', params);
Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params);

params.C_h_real = 1/2*eye(2*M*K);
params.Phi = Phi;

nmse_blmmse = zeros(size(SNRs));
nmse_blmmse_cn = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_blmmse = 0;
    current_error_blmmse_cn = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    params_blmmse = mtk_prepare_channel_estimation('blmmse', params);
    params_blmmse_cn = mtk_prepare_channel_estimation('blmmse-cn', params);

    for channel_index=1:iterations
        params_blmmse.H = H(:,:,channel_index);
        params_blmmse_cn.H = H(:,:,channel_index);
        params_blmmse.N = N(:,:,1);
        params_blmmse_cn.N = N(:,:,1);

        H_hat = mtk_estimate_channel('blmmse', params_blmmse);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse = current_error_blmmse + norm(E, 'fro')^2;

        H_hat = mtk_estimate_channel('blmmse-cn', params_blmmse_cn);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse_cn = current_error_blmmse_cn + norm(E, 'fro')^2;
    end
    nmse_blmmse(snr_index) = current_error_blmmse/(iterations*M*K);
    nmse_blmmse_cn(snr_index) = current_error_blmmse_cn/(iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

plot(SNRs, nmse_blmmse, 'r-o', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, nmse_blmmse_cn, 'b-o', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

legend('Numerical MSE_{BLMMSE}', 'Numerical MSE_{BLMMSE-CN}');
xlabel('\rho [dB]')
ylim([0.3, 0.75]);
ylabel('MSE of channel estimation')
