params.K = 3;
params.M = 16;
params.tau = 3;
params.channel_iterations = 1000;
params.channel_index = 11;
params.noise_iterations = params.channel_iterations * params.channel_index;
params.SNR = -5;
params.r_k = 0.5;
params.eta = mtk_util_jakes(3/3.6, 2.5e9, 5e-3)*ones(params.K, 1);

tic;
[params.H, results] = mtk_generate_channel('kron-markov', params, 0);
params.Phi = mtk_generate_pilot('dft', params);
params.N = mtk_generate_noise(params, 3);

params.C_h_real = results.C_h_real;
params.eta_real = results.eta_real;
params.zeta_real = results.zeta_real;

sim.quantizer = '1bit';
sim.params = params;

tic;
sim.method = 'kfb';
[nmse, results] = mtk_sim_ce_index(sim);

params.alpha = 32;
[~, params.Beff] = mtk_util_comparator_network('random', params, 0);
sim.method = 'kfb-cn';
sim.params = params;
[nmse_cn, results_cn] = mtk_sim_ce_index(sim);
toc;

plot(0:params.channel_index-1, 10*log10(nmse), 'b-o', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(0:params.channel_index-1, 10*log10(results.mmse), 'b--', 'LineWidth', 1, 'MarkerSize', 8);

plot(0:params.channel_index-1, 10*log10(nmse_cn), 'r-o', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(0:params.channel_index-1, 10*log10(results_cn.mmse), 'r--', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

xlabel('Time Index')
ylabel('NMSE [dB]')
xlim([0 params.channel_index-1]);
ylim([-7.5 -1.5]);
grid();
