tic;
channel_iterations = 1000;
K = 2;
M = 2;
tau = 16;
alpha = 4;

SNRs = linspace(-10, 20, 21);

params.K = K;
params.M = M;
params.tau = tau;
params.channel_iterations = channel_iterations;
params.noise_iterations = channel_iterations;
params.alpha = alpha;

[params.B, params.Beff] = mtk_util_comparator_network('random', params, 0);

[H, results] = mtk_generate_channel('rayleigh', params, 1);
params.C_h_real = results.C_h_real;
params.Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params, 2);

params.quantizer = '1bit';

nmse_bmmse = zeros(size(SNRs));
nmse_bmmse_cn = zeros(size(SNRs));
nmse_ls = zeros(size(SNRs));
nmse_aqn = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_bmmse_cn = 0;
    current_error_bmmse = 0;
    current_error_ls = 0;
    current_error_aqn = 0;

    params.rho = mtk_util_db_to_linear(SNR);

    params_bmmse_cn = mtk_ce_prepare('bmmse-cn', params);
    params_bmmse = mtk_ce_prepare('bmmse', params);
    params_ls = mtk_ce_prepare('ls', params);
    params_aqn = mtk_ce_prepare('aqn', params);

    for channel_index=1:channel_iterations
        params_bmmse_cn.H = H(:,:,channel_index);
        params_bmmse.H = H(:,:,channel_index);
        params_ls.H = H(:,:,channel_index);
        params_aqn.H = H(:,:,channel_index);
        
        params_bmmse_cn.N = N(:,:,channel_index);
        params_bmmse.N = N(:,:,channel_index);
        params_ls.N = N(:,:,channel_index);
        params_aqn.N = N(:,:,channel_index);

        H_hat = mtk_ce_estimate('bmmse-cn', params_bmmse_cn);
        E = H(:,:,channel_index) - H_hat;
        current_error_bmmse_cn = current_error_bmmse_cn + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('bmmse', params_bmmse);
        E = H(:,:,channel_index) - H_hat;
        current_error_bmmse = current_error_bmmse + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('ls', params_ls);
        E = H(:,:,channel_index) - H_hat;
        current_error_ls = current_error_ls + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('aqn', params_aqn);
        E = H(:,:,channel_index) - H_hat;
        current_error_aqn = current_error_aqn + norm(E, 'fro')^2;
    end
    nmse_bmmse_cn(snr_index) = current_error_bmmse_cn/(channel_iterations*M*K);
    nmse_bmmse(snr_index) = current_error_bmmse/(channel_iterations*M*K);
    nmse_ls(snr_index) = current_error_ls/(channel_iterations*M*K);
    nmse_aqn(snr_index) = current_error_aqn/(channel_iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

figure(1);
plot(SNRs, nmse_ls, '-->');
hold on;
plot(SNRs, nmse_bmmse, 'r-s');
hold on;
plot(SNRs, nmse_bmmse_cn, 'b-s');
plot(SNRs, nmse_aqn, '--x');
hold off;

legend('LS', 'BMMSE', 'BMMSE-CN', 'Additive Quantizer Noise');
xlabel('SNR [dB]')
ylabel('NMSE')
