function x = mtk_generate_symbols(modulation, params, seed)
    switch modulation
        case 'qpsk'
            K = params.K;
            symbol_iterations = params.symbol_iterations;

            alphabet = 1/sqrt(2)*[1+1i -1+1i -1-1i 1-1i];

            rng(seed);
            x = alphabet(randi(length(alphabet), K, symbol_iterations));
        otherwise
            error('Modulation not implemented');
    end
end

