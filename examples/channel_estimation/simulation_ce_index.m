clear all;

params.K = 2;
params.M = 2;
params.tau = 2;
params.channel_iterations = 10;
params.noise_iterations = 30;
params.channel_index = 30;
params.SNR = -5;
params.r_k = 0.8;
params.eta = 0.998*ones(params.K, 1);

tic;
params.H = mtk_generate_channel('kron-markov', params, 0);
params.Phi = mtk_generate_pilot('dft', params);
params.N = mtk_generate_noise(params, 1);

% TODO: get from generation
params.C_h_real = 1/2*eye(2 * params.M * params.K);

sim.quantizer = '1bit';
sim.params = params;

tic;
sim.method = 'kfb';
nmse = mtk_sim_ce_index(sim);
toc;

plot(1:params.channel_index, 10*log10(nmse), '-->', 'Color', '#da7e26', 'LineWidth', 1, 'MarkerSize', 8);
hold on;

xlabel('Time Index')
ylabel('NMSE [dB]')
