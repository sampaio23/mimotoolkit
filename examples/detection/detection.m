tic;
channel_iterations = 100;
symbol_iterations = 100;
noise_iterations = 100;
K = 2;
M = 2;
tau = 1;

SNRs = linspace(-20, 20, 30);

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

error_zf = zeros(size(SNRs));
error_mrc = zeros(size(SNRs));
error_mmse = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_zf = 0;
    current_error_mrc = 0;
    current_error_mmse = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    for channel_index=1:channel_iterations
        params.H = H(:,:,channel_index);
        params.H_hat = params.H; % Perfect CSI
        G_zf = mtk_det_detector('zf', params);
        G_mrc = mtk_det_detector('mrc', params);
        G_mmse = mtk_det_detector('mmse', params);
        for noise_index=1:noise_iterations
            params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);

            r = sqrt(params.rho) * params.H * params.x + params.n;
            y = r;

            x_tilde = G_zf * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_zf = current_error_zf + sum(x_hat ~= params.x, 'all');

            x_tilde = G_mrc * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_mrc = current_error_mrc + sum(x_hat ~= params.x, 'all');

            x_tilde = G_mmse * y;
            x_hat = mtk_util_demodulate('qpsk', x_tilde);
            current_error_mmse = current_error_mmse + sum(x_hat ~= params.x, 'all');
        end
    end
    error_zf(snr_index) = current_error_zf / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_mrc(snr_index) = current_error_mrc / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    error_mmse(snr_index) = current_error_mmse / (params.K * channel_iterations * symbol_iterations * noise_iterations);
    snr_index = snr_index + 1;
end
toc;

semilogy(SNRs, error_zf, '-->');
hold on;
semilogy(SNRs, error_mrc, '--o');
semilogy(SNRs, error_mmse, '--o');
hold off;

legend("ZF", "MRC", "MMSE");
xlim([min(SNRs) max(SNRs)]);

%%
scatter(real(x(1,:)), imag(x(1,:)));
hold on;
scatter(real(x(2,:)), imag(x(2,:)));
hold off;