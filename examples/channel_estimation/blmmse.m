clear all
tic;
channel_iterations = 10;
noise_iterations = 100;
K = 2;
M = 2;
tau = 2;
alpha = 2;

SNRs = linspace(-10, 20, 21);

params.K = K;
params.M = M;
params.tau = tau;
params.channel_iterations = channel_iterations;
params.noise_iterations = noise_iterations;
params.alpha = alpha;

[params.B, params.Beff] = mtk_util_comparator_network('random', params, 0);

params.r_k = 0.5;
params.channel_index = 1;
params.eta = rand(1,K);
[H, results] = mtk_generate_channel('kron-markov', params, 0);
Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params, 1);

params.C_h_real = results.C_h_real;
params.Phi = Phi;
params.quantizer = '1bit';

nmse_bmmse = zeros(size(SNRs));
nmse_bmmse_cn = zeros(size(SNRs));
nmse_ls = zeros(size(SNRs));
nmse_adn = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_bmmse_cn = 0;
    current_error_bmmse = 0;
    current_error_ls = 0;
    current_error_adn = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    params_bmmse_cn = mtk_ce_prepare('blmmse-cn', params);
    params_bmmse = mtk_ce_prepare('blmmse', params);
    params_ls = mtk_ce_prepare('ls', params);
    params_adn = mtk_ce_prepare('aqnm', params);

    for channel_index=1:channel_iterations
        params_bmmse_cn.H = H(:,:,channel_index);
        params_bmmse.H = H(:,:,channel_index);
        params_ls.H = H(:,:,channel_index);
        params_adn.H = H(:,:,channel_index);
        for noise_index=1:noise_iterations
            params_bmmse_cn.N = N(:,:,noise_index);
            params_bmmse.N = N(:,:,noise_index);
            params_ls.N = N(:,:,noise_index);
            params_adn.N = N(:,:,noise_index);
    
            H_hat = mtk_ce_estimate('blmmse-cn', params_bmmse_cn);
            E = H(:,:,channel_index) - H_hat;
            current_error_bmmse_cn = current_error_bmmse_cn + norm(E, 'fro')^2;
    
            H_hat = mtk_ce_estimate('blmmse', params_bmmse);
            E = H(:,:,channel_index) - H_hat;
            current_error_bmmse = current_error_bmmse + norm(E, 'fro')^2;
    
            H_hat = mtk_ce_estimate('ls', params_ls);
            E = H(:,:,channel_index) - H_hat;
            current_error_ls = current_error_ls + norm(E, 'fro')^2;
    
            H_hat = mtk_ce_estimate('aqnm', params_adn);
            E = H(:,:,channel_index) - H_hat;
            current_error_adn = current_error_adn + norm(E, 'fro')^2;
        end
    end
    nmse_bmmse_cn(snr_index) = current_error_bmmse_cn/(noise_iterations*channel_iterations*M*K);
    nmse_bmmse(snr_index) = current_error_bmmse/(noise_iterations*channel_iterations*M*K);
    nmse_ls(snr_index) = current_error_ls/(noise_iterations*channel_iterations*M*K);
    nmse_adn(snr_index) = current_error_adn/(noise_iterations*channel_iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

plot(SNRs, nmse_ls, '-->', 'Color', '#da7e26', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, nmse_bmmse, 'r-s', 'Color', '#eb1f24', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, nmse_bmmse_cn, 'b-s', 'Color', '#0000ff', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, nmse_adn, '--x', 'Color', '#017f3f', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

legend('LS','BMMSE', 'BMMSE-CN', 'Additive Quantizer Noise');
xlabel('SNR (dB)')
ylabel('Normalized MSE (dB)')
ylabel('MSE')
