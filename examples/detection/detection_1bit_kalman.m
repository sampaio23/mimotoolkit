tic;
channel_iterations = 1000;
symbol_iterations = 100;
channel_index = 31;
noise_iterations = 100;
K = 2;
M = 16;
tau = 2;

SNRs = linspace(-10, 30, 1);

params.K = K;
params.M = M;
params.channel_iterations = channel_iterations;
params.symbol_iterations = symbol_iterations;
params.channel_index = channel_index;
params.noise_iterations = noise_iterations;
params.quantizer = '1bit';
params.r_k = 0.5;
params.eta = mtk_util_jakes(3/3.6, 2.5e9, 5e-3)*ones(params.K, 1);

[H, results] = mtk_generate_channel('kron-markov', params, 0);
params.C_h_real = results.C_h_real;
params.eta_real = results.eta_real;
params.zeta_real = results.zeta_real;
params.x = mtk_generate_symbols('qpsk', params, 1);

params.tau = 1;
n = mtk_generate_noise(params, 2);

params.tau = tau;
params.Phi = mtk_generate_pilot('dft', params);
params.noise_iterations = channel_iterations*channel_index;
N = mtk_generate_noise(params, 3);

error_bmmse = zeros(size(SNRs));

snr_index = 1;
for SNR = SNRs
    current_error_bmmse = 0;
    params.rho = mtk_util_db_to_linear(SNR);
    params = mtk_ce_prepare('kfb', params);
    for it=1:channel_iterations
        params.H = H(:,:,:,it);
        params.N = N(:,:,(it-1)*channel_index+1:it*channel_index);

        H_hat = mtk_ce_estimate('kfb', params); % Imperfect CSI
        params.H_hat = H_hat(:,:,channel_index);
        G_zf = mtk_det_detector('bmmse', params);
        params.H = H(:,:,channel_index,it);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);
    
            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = mtk_util_quantize('1bit', r, 1);

            x_tilde = G_zf * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');
        end
    end
    error_bmmse(snr_index) = current_error_bmmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end

error_bmmse_rob = zeros(size(SNRs));

snr_index = 1;
for SNR = SNRs
    current_error_bmmse = 0;
    params.rho = mtk_util_db_to_linear(SNR);
    params = mtk_ce_prepare('kfb', params);
    for it=1:channel_iterations
        params.H = H(:,:,:,it);
        params.N = N(:,:,(it-1)*channel_index+1:it*channel_index);

        [H_hat, results] = mtk_ce_estimate('kfb', params); % Imperfect CSI
        params.H_hat = H_hat(:,:,channel_index);
        params.error_real = results.mmse(:,:,channel_index);
        Gamma_real = mtk_util_ce_gamma('qpsk', params);
        params.Gamma = mtk_util_mat_complex(Gamma_real);
        G_zf = mtk_det_detector('robust-bmmse', params);
        params.H = H(:,:,channel_index,it);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);
    
            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = mtk_util_quantize('1bit', r, 1);

            x_tilde = G_zf * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');
        end
    end
    error_bmmse_rob(snr_index) = current_error_bmmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end
toc;

semilogy(SNRs, error_bmmse/2, '-o');
hold on;
semilogy(SNRs, error_bmmse_rob/2, '-o');

legend("BMMSE", "Robust BMMSE");
xlim([min(SNRs) max(SNRs)]);
ylim([1e-3 5e-1]);
grid();