function x = mtk_generate_symbols(modulation, params)
    switch modulation
        case 'qpsk'
            K = params.K;
            symbol_iterations = params.symbol_iterations;

            alphabet = 1/sqrt(2)*[1+1i -1+1i -1-1i 1-1i];
            
            current_seed = rng;
            rng("shuffle");
            x = alphabet(randi(length(alphabet), K, symbol_iterations));
            rng(current_seed);
        otherwise
            error('Modulation not implemented');
    end
end

