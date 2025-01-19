channel_iterations = 10;
noise_iterations = 1;
K = 4;
M = 16;
tau = 20;

SNRs = linspace(-20, 20, 21);

params.K = K;
params.M = M;
params.tau = tau;
params.channel_iterations = channel_iterations;
params.noise_iterations = 10;

params.H = mtk_generate_channel('rayleigh', params, 0);
params.Phi = mtk_generate_pilot('dft', params);
params.N = mtk_generate_noise(params, 1);

params.C_h_real = 1/2*eye(2*M*K);

sim.quantizer = '1bit';
sim.SNRs = SNRs;
sim.params = params;

tic;
sim.method = 'ls';
nmse_ls = mtk_sim_ce_snr(sim);
sim.method = 'blmmse';
nmse_blmmse = mtk_sim_ce_snr(sim);
sim.method = 'aqnm';
nmse_aqnm = mtk_sim_ce_snr(sim);
toc;

plot(SNRs, 10*log10(nmse_ls), '-->', 'Color', '#da7e26', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, 10*log10(nmse_blmmse), 'r-s', 'Color', '#eb1f24', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(SNRs, 10*log10(nmse_aqnm), '--x', 'Color', '#017f3f', 'LineWidth', 1, 'MarkerSize', 8);
hold off;

legend('LS Proposed in [21]','BLMMSE', 'Additive Quantizer Noise in [28]');
xlabel('SNR (dB)')
ylim([-8,8]);
ylabel('Normalized MSE (dB)')
