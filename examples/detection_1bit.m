tic;
channel_iterations = 100;
symbol_iterations = 100;
noise_iterations = 100;
K = 2;
M = 5;
tau = 1;

SNRs = linspace(-5, 30, 8);

params.K = K;
params.M = M;
params.iterations = channel_iterations;
params.symbol_iterations = symbol_iterations;
params.noise_iterations = noise_iterations;
params.tau = tau;

H = mtk_generate_channel('rayleigh', params);
x = mtk_generate_symbols('qpsk', params);
n = mtk_generate_noise(params);

params.x = x;

error_mrc = zeros(size(SNRs));
error_bmrc = zeros(size(SNRs));
error_zf = zeros(size(SNRs));
error_bzf = zeros(size(SNRs));
error_mmse = zeros(size(SNRs));
error_bmmse = zeros(size(SNRs));

snr_index = 1;
for SNR = SNRs
    current_error_mrc = 0;
    current_error_bmrc = 0;
    current_error_zf = 0;
    current_error_bzf = 0;
    current_error_mmse = 0;
    current_error_bmmse = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    for channel_index=1:channel_iterations
        params.H = H(:,:,channel_index);
        params.H_hat = params.H; % Perfect CSI

        G_mrc = mtk_detector('mrc', params);
        G_bmrc = mtk_detector('bmrc', params);
        G_zf = mtk_detector('zf', params);
        G_bzf = mtk_detector('bzf', params);
        G_mmse = mtk_detector('mmse', params);
        G_bmmse = mtk_detector('bmmse', params);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);
    
            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = mtk_util_quantize(r, (params.K*params.rho + 1)/(pi*params.rho));

            x_tilde = G_mrc * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_mrc = current_error_mrc + sum(x_hat ~= params.x, 'all');

            x_tilde = G_bmrc * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmrc = current_error_bmrc + sum(x_hat ~= params.x, 'all');

            x_tilde = G_zf * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_zf = current_error_zf + sum(x_hat ~= params.x, 'all');

            x_tilde = G_bzf * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bzf = current_error_bzf + sum(x_hat ~= params.x, 'all');

            x_tilde = G_mmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_mmse = current_error_mmse + sum(x_hat ~= params.x, 'all');

            x_tilde = G_bmmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');
        end
    end
    error_mrc(snr_index) = current_error_mrc / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bmrc(snr_index) = current_error_bmrc / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_zf(snr_index) = current_error_zf / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bzf(snr_index) = current_error_bzf / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_mmse(snr_index) = current_error_mmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_bmmse(snr_index) = current_error_bmmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end
toc;

semilogy(SNRs, error_mrc/2, 'b--');
hold on;
semilogy(SNRs, error_bmrc/2, 'r-');
semilogy(SNRs, error_zf/2, 'b--*');
semilogy(SNRs, error_bzf/2, 'r-*');
semilogy(SNRs, error_mmse/2, 'b--o');
semilogy(SNRs, error_bmmse/2, 'r-o');
hold off;

legend("MRC", "BMRC", "ZF", "BZF", "MMSE", "BMMSE");
xlim([min(SNRs) max(SNRs)]);

%%
scatter(real(x(1,:)), imag(x(1,:)));
hold on;
scatter(real(x(2,:)), imag(x(2,:)));
hold off;