tic;
channel_iterations = 10;
symbol_iterations = 11;
noise_iterations = 9;
K = 2;
M = 16;
tau = 1;

SNRs = linspace(-5, 15, 8);

params.K = K;
params.M = M;
params.channel_iterations = channel_iterations;
params.symbol_iterations = symbol_iterations;
params.noise_iterations = noise_iterations;
params.tau = tau;

params.H = mtk_generate_channel('rayleigh', params, 0);
params.x = mtk_generate_symbols('qpsk', params, 1);
params.n = mtk_generate_noise(params, 2);

sim.method = 'perfect';
sim.modulation = 'qpsk';
sim.detector = 'bmmse';
sim.quantizer = '1bit';
sim.SNRs = SNRs;
sim.params = params;

SER = mtk_sim_detect_snr(sim);
toc;

figure();
semilogy(SNRs, SER, 'r-o');
xlim([min(SNRs) max(SNRs)]);
