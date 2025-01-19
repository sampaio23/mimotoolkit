function H = mtk_generate_channel(model, params, seed)
    switch model
        case 'rayleigh'
            K = params.K;
            M = params.M;
            iterations = params.channel_iterations;

            rng(seed);
            H = 1/sqrt(2)*(randn(M, K, iterations) + 1i*randn(M, K, iterations));
        otherwise
            error('Channel type not implemented');
    end
end
