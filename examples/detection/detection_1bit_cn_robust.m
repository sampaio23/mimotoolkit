tic;
channel_iterations = 1000;
symbol_iterations = 100;
noise_iterations = 100;
K = 4;
M = 16;
tau = 4;

SNRs = linspace(-20, 20, 9);

params.K = K;
params.M = M;
params.channel_iterations = channel_iterations;
params.symbol_iterations = symbol_iterations;
params.noise_iterations = noise_iterations;
params.quantizer = '1bit';

[H, results] = mtk_generate_channel('rayleigh', params, 0);
params.C_h_real = results.C_h_real;
params.x = mtk_generate_symbols('qpsk', params, 1);

params.tau = 1;
n = mtk_generate_noise(params, 2);

params.tau = tau;
params.Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params, 3);

params.alpha = 2*M;
[params.B, params.Beff] = mtk_util_comparator_network('random', params, 4);

error_bmmse = zeros(size(SNRs));
error_bmmse_imp = zeros(size(SNRs));
error_bmmse_rob = zeros(size(SNRs));

snr_index = 1;
for SNR = SNRs
    current_error_bmmse = 0;
    current_error_bmmse_imp = 0;
    current_error_bmmse_rob = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    [params, results] = mtk_ce_prepare('blmmse', params);
    params.error_real = results.error_real;
    Gamma_real = mtk_util_ce_gamma('qpsk', params);
    params.Gamma = mtk_util_mat_complex(Gamma_real); % Can I do this? Test with CN
    for channel_index=1:channel_iterations
        params.H = H(:,:,channel_index);
        params.N = N(:,:,1);
        H_hat = params.H; % Perfect CSI
        H_hat_imp = mtk_ce_estimate("blmmse", params); % Imperfect CSI

        params.H_hat = H_hat;
        G_bmmse = mtk_det_detector('bmmse', params);
        params.H_hat = H_hat_imp;
        G_bmmse_imp = mtk_det_detector('bmmse', params);
        G_bmmse_rob = mtk_det_detector('robust-bmmse', params);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);
    
            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = mtk_util_quantize('1bit', r, 1);

            x_tilde = G_bmmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');

            x_tilde = G_bmmse_imp * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse_imp = current_error_bmmse_imp + sum(x_hat ~= params.x, 'all');

            x_tilde = G_bmmse_rob * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse_rob = current_error_bmmse_rob + sum(x_hat ~= params.x, 'all');
        end
    end
    error_bmmse(snr_index) = current_error_bmmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bmmse_imp(snr_index) = current_error_bmmse_imp / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bmmse_rob(snr_index) = current_error_bmmse_rob / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end
toc;

semilogy(SNRs, error_bmmse/2, 'r--o');
hold on;
semilogy(SNRs, error_bmmse_imp/2, 'r-o');
semilogy(SNRs, error_bmmse_rob/2, 'b-o');
hold off;

legend("Perfect CSI", "Imperfect CSI", "Robust");
xlim([min(SNRs) max(SNRs)]);
ylim([5e-4 1e0]);
grid();