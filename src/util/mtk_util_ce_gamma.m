function Gamma_real = mtk_util_ce_gamma(modulation, params)
    switch modulation
        case 'qpsk'
            alphabet = 1/sqrt(2)*[1+1i -1+1i -1-1i 1-1i]; % TODO: create 
            % get alphabet to avoid copy in mtk_generate_symbols
        otherwise
            error('Modulation not implemented');
    end

    X = cell(1, params.K);
    [X{:}] = ndgrid(alphabet);
    X = reshape(cat(params.K+1, X{:}), [], params.K);
    
    Gamma_real = zeros(2*params.M, 2*params.M);
    for i=1:length(alphabet)^params.K
        X_tilde = kron(X(i,:), eye(params.M));
        X_tilde_real = mtk_util_mat_real(X_tilde);
        Gamma_real = Gamma_real + X_tilde_real * params.error_real * X_tilde_real';
    end
    Gamma_real = Gamma_real / length(alphabet)^params.K;
end

