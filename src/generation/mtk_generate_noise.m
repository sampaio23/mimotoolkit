function N = mtk_generate_noise(params, seed)
    tau = params.tau;
    M = params.M;
    noise_iterations = params.noise_iterations;

    rng(seed);
    N = 1/sqrt(2)*(randn(M, tau, noise_iterations) + 1i*randn(M, tau, noise_iterations));
end

