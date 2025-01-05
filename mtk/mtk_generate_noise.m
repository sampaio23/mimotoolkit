function N = mtk_generate_noise(params)
    tau = params.tau;
    M = params.M;
    iterations = params.iterations;

    current_seed = rng(1009*iterations);
    N = 1/sqrt(2)*(randn(M, tau, iterations) + 1i*randn(M, tau, iterations));
    rng(current_seed);
end

