tic;
channel_iterations = 10;
noise_iterations = 1;
K = 2;
M = 2;
tau = 2;
alpha = 1;

SNRs = linspace(-10, 20, 9);

params.K = K;
params.M = M;
params.tau = tau;
params.channel_iterations = channel_iterations;
params.noise_iterations = noise_iterations;

params.alpha = alpha;
[params.B, params.Beff] = mtk_util_comparator_network('random', params, 0);

H = mtk_generate_channel('rayleigh', params, 0);
Phi = mtk_generate_pilot('dft', params);
N = mtk_generate_noise(params, 1);

params.C_h_real = 1/2*eye(2*M*K);
params.Phi = Phi;
params.quantizer = '1bit';

nmse_blmmse = zeros(size(SNRs));
nmse_blmmse_analytical = zeros(size(SNRs));
nmse_blmmse_cn = zeros(size(SNRs));
nmse_blmmse_cn_analytical = zeros(size(SNRs));
snr_index = 1;
for SNR = SNRs
    current_error_blmmse = 0;
    current_error_blmmse_cn = 0;

    params.rho = mtk_util_db_to_linear(SNR);
    params_blmmse = mtk_ce_prepare('blmmse', params);
    params_blmmse_cn = mtk_ce_prepare('blmmse-cn', params);

    nmse_blmmse_analytical(snr_index) = mtk_ce_estimate_analytical('blmmse', params_blmmse);
    nmse_blmmse_cn_analytical(snr_index) = mtk_ce_estimate_analytical('blmmse-cn', params_blmmse_cn);

    for channel_index=1:channel_iterations
        params_blmmse.H = H(:,:,channel_index);
        params_blmmse_cn.H = H(:,:,channel_index);
        params_blmmse.N = N(:,:,1);
        params_blmmse_cn.N = N(:,:,1);

        H_hat = mtk_ce_estimate('blmmse', params_blmmse);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse = current_error_blmmse + norm(E, 'fro')^2;

        H_hat = mtk_ce_estimate('blmmse-cn', params_blmmse_cn);
        E = H(:,:,channel_index) - H_hat;
        current_error_blmmse_cn = current_error_blmmse_cn + norm(E, 'fro')^2;
    end
    nmse_blmmse(snr_index) = current_error_blmmse/(channel_iterations*M*K);
    nmse_blmmse_cn(snr_index) = current_error_blmmse_cn/(channel_iterations*M*K);
    snr_index = snr_index + 1;
end
toc;

plot(SNRs, nmse_blmmse_analytical, 'r-', 'LineWidth', 1);
hold on;
plot(SNRs, nmse_blmmse, 'ro', 'LineWidth', 1, 'MarkerSize', 8);
plot(SNRs, nmse_blmmse_cn_analytical, 'b-', 'LineWidth', 1);
plot(SNRs, nmse_blmmse_cn, 'bo', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

legend('Analytical MSE_{BLMMSE}', 'Numerical MSE_{BLMMSE}', 'Analytical MSE_{BLMMSE-CN}', 'Numerical MSE_{BLMMSE-CN}');
xlabel('\rho [dB]')
% ylim([0.3, 0.75]);
ylabel('MSE of channel estimation')

%%
SNRs = linspace(-10, 20, 9);

nmse_blmmse_cn_nn = [2.518674, 2.421271, 2.160343, 1.902718, 1.636279, 1.511470, 1.442313, 1.418068, 1.401485];
nmse_blmmse_cn_nn = nmse_blmmse_cn_nn/(M*K);
nmse_blmmse_cn_nn_2 = [2.517816, 2.403144, 2.157796, 1.867760, 1.608106, 1.464774, 1.388171, 1.368336, 1.355259];
nmse_blmmse_cn_nn_2 = nmse_blmmse_cn_nn_2/(M*K);

plot(SNRs, nmse_blmmse_analytical, 'r-o', 'LineWidth', 1);
hold on;
plot(SNRs, nmse_blmmse_cn_analytical, 'b-o', 'LineWidth', 1);
plot(SNRs, nmse_blmmse_cn_nn, 'k-o', 'LineWidth', 1);
plot(SNRs, nmse_blmmse_cn_nn_2, 'g-o', 'LineWidth', 1);
hold off;

legend('NMSE_{BLMMSE}', 'NMSE_{BLMMSE-CN} 1', 'NMSE_{BLMMSE-CN-NN} 1', 'NMSE_{BLMMSE-CN-NN} 2');
xlabel('\rho [dB]')
% ylim([0.3, 0.75]);
ylabel('MSE of channel estimation')

%%
SNRs = linspace(-10, 20, 9);

nmse_blmmse_cn_nn = [2.403562 2.160607 1.903442 1.615632 1.419671 1.312652 1.264963 1.245861 1.231686];
nmse_blmmse_cn_nn = nmse_blmmse_cn_nn/(M*K);

plot(SNRs, nmse_blmmse_analytical, 'r-o', 'LineWidth', 1);
hold on;
plot(SNRs, nmse_blmmse_cn_analytical, 'b-o', 'LineWidth', 1);
plot(SNRs, nmse_blmmse_cn_nn, 'k-o', 'LineWidth', 1);
hold off;

legend('Numerical BLMMSE', 'Numerical BLMMSE-CN Independent', 'Neural');
xlabel('\rho [dB]')
% ylim([0.3, 0.75]);
ylabel('MSE of channel estimation')