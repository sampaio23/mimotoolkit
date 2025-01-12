function N = mtk_generate_noise(params)
    tau = params.tau;
    M = params.M;
    noise_iterations = params.noise_iterations;

    current_seed = rng;
    rng("shuffle");
    N = 1/sqrt(2)*(randn(M, tau, noise_iterations) + 1i*randn(M, tau, noise_iterations));
    rng(current_seed);
end

