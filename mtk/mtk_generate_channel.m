function H = mtk_generate_channel(model, params)
    switch model
        case 'rayleigh'
            K = params.K;
            M = params.M;
            iterations = params.iterations;

            current_seed = rng(2003*iterations);
            H = 1/sqrt(2)*(randn(M, K, iterations) + 1i*randn(M, K, iterations));
            rng(current_seed);
        otherwise
            error('Channel type not implemented');
    end
end
