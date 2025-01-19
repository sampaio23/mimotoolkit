function H_hat = mtk_ce_estimate(method, params)
    switch method
        case 'perfect'
            H_hat = params.H;
        case 'ls'
            H = params.H;
            Phi = params.Phi;
            N = params.N;
            rho = params.rho;

            Y = sqrt(rho)*H*conj(Phi') + N;
            R = mtk_util_quantize(params.quantizer, Y, 1/sqrt(2));

            H_hat = (1/sqrt(rho))*(R*conj(Phi))/(conj(Phi')*conj(Phi));
        case {'blmmse','aqnm'}
            C_h_real = params.C_h_real;
            H = params.H;
            K = params.K;
            M = params.M;
            N = params.N;
            Phi = params.Phi;
            rho = params.rho;

            Phi_tilde_real = params.Phi_tilde_real;
            C_r_real = params.C_r_real;

            Y = sqrt(rho)*H*conj(Phi') + N;
            
            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, y_real, 1);
            
            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat = reshape(h_hat, [M, K]);
        case 'blmmse-cn'
            Beff = params.Beff;
            C_h_real = params.C_h_real;
            H = params.H;
            K = params.K;
            M = params.M;
            N = params.N;
            Phi = params.Phi;
            rho = params.rho;

            Phi_tilde_real = params.Phi_tilde_real;
            C_r_real = params.C_r_real;

            Y = sqrt(rho)*H*conj(Phi') + N;

            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, Beff * y_real, 1);

            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat = reshape(h_hat, [M, K]);
        otherwise
            error('Channel type not implemented');
    end
end

