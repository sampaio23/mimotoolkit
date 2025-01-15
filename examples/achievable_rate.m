tic;
channel_iterations = 1000;
K = 3;
M = 8;

SNRs = linspace(-20, 30, 11);

params.K = K;
params.M = M;
params.tau = K;
params.iterations = channel_iterations;

H = mtk_generate_channel('rayleigh', params);

params.alpha = M;
[params.B, ~] = mtk_util_comparator_network('independent', params);

rate = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_rate = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    for channel_index=1:channel_iterations
        params.H = H(:,:,channel_index);
        params.H_hat = params.H; % Perfect CSI

        I = mtk_achievable_rate(params);
        current_rate = current_rate + I;
    end
    rate(snr_index) = current_rate / (channel_iterations);
    snr_index = snr_index + 1;
end

toc;

figure();
plot(SNRs, rate, 'r-o');
grid();
