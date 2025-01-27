tic;
channel_iterations = 1000;
symbol_iterations = 100;
noise_iterations = 100;
K = 4;
M = 16;
tau = 1;

SNRs = linspace(-10, 30, 9);

params.K = K;
params.M = M;
params.channel_iterations = channel_iterations;
params.symbol_iterations = symbol_iterations;
params.noise_iterations = noise_iterations;
params.tau = tau;

H = mtk_generate_channel('rayleigh', params, 0);
x = mtk_generate_symbols('qpsk', params, 1);
n = mtk_generate_noise(params, 2);

params.x = x;

params.alpha = 2*M;
[params.B, params.Beff] = mtk_util_comparator_network('random', params, 4);

error_bmmse = zeros(size(SNRs));
error_bmmse_cn = zeros(size(SNRs));
error_mmse_unquant = zeros(size(SNRs));

snr_index = 1;
for SNR = SNRs
    current_error_bmmse = 0;
    current_error_bmmse_cn = 0;
    current_error_mmse_unquant = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    for channel_index=1:channel_iterations
        params.H = H(:,:,channel_index);
        params.H_hat = params.H; % Perfect CSI

        params.x_real = mtk_util_vec_real(params.x);

        G_bmmse = mtk_det_detector('bmmse', params);
        G_mmse = mtk_det_detector('mmse', params);
        G_bmmse_cn = mtk_det_detector('bmmse-cn', params);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);
    
            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = mtk_util_quantize('1bit', r, 1);
            x_tilde = G_bmmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');

            y = mtk_util_quantize('unquant', r, 1);
            x_tilde = G_mmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_mmse_unquant = current_error_mmse_unquant + sum(x_hat ~= params.x, 'all');

            params.H_real = mtk_util_mat_real(params.H);
            params.n_real = mtk_util_vec_real(params.n);
            r_real = params.B * (sqrt(params.rho) * params.H_real * params.x_real + params.n_real);
            y_real = mtk_util_quantize('1bit', r_real, 1);
            x_tilde_real = G_bmmse_cn * y_real;
            x_tilde = mtk_util_vec_complex(x_tilde_real);
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse_cn = current_error_bmmse_cn + sum(x_hat ~= params.x, 'all');
        end
    end
    error_bmmse(snr_index) = current_error_bmmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bmmse_cn(snr_index) = current_error_bmmse_cn / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_mmse_unquant(snr_index) = current_error_mmse_unquant / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end
toc;

semilogy(SNRs, error_bmmse/2, 'r-o');
hold on;
semilogy(SNRs, error_bmmse_cn/2, 'b-o');
semilogy(SNRs, error_mmse_unquant/2, 'k-*');
hold off;
xlim([min(SNRs) max(SNRs)]);
ylim([1e-5 5e-1]);
grid();
