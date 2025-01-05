tic;
iterations = 1000;
K = 4;
M = 16;
tau = 20;

SNRs = linspace(-20, 20, 21);

params.K = K;
params.M = M;
params.tau = tau;
params.iterations = iterations;

H = mtk_generate_channel('rayleigh', params);
Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params);

nmse_blmsse = zeros(size(SNRs));
nmse_ls = zeros(size(SNRs));
nmse_adn = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_blmmse = 0;
    current_error_ls = 0;
    current_error_adn = 0;
    for channel_index=1:iterations
        params.H = H(:,:,channel_index);
        params.Phi = Phi;
        params.N = N(:,:,1);
        params.rho = mtk_util_db_to_linear(SNR);
        params.C_h_real = 1/2*eye(2*M*K);

        H_hat = mtk_estimate_channel('blmmse', params);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse = current_error_blmmse + norm(E, 'fro')^2;

        H_hat = mtk_estimate_channel('ls', params);
        E = H(:,:,channel_index) - H_hat;
        current_error_ls = current_error_ls + norm(E, 'fro')^2;

        H_hat = mtk_estimate_channel('additive-quantizer-noise', params);
        E = H(:,:,channel_index) - H_hat;
        current_error_adn = current_error_adn + norm(E, 'fro')^2;
    end
    nmse_blmsse(snr_index) = current_error_blmmse/(iterations*M*K);
    nmse_ls(snr_index) = current_error_ls/(iterations*M*K);
    nmse_adn(snr_index) = current_error_adn/(iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

%%
plot(SNRs, 10*log10(nmse_ls), '-->', 'Color', '#da7e26', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, 10*log10(nmse_blmsse), 'r-s', 'Color', '#eb1f24', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, 10*log10(nmse_adn), '--x', 'Color', '#017f3f', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

legend('LS Proposed in [21]','BLMMSE', 'Additive Quantizer Noise in [28]');
xlabel('SNR (dB)')
ylim([-8,8]);
ylabel('Normalized MSE (dB)')
hold on;
