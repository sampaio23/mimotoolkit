function H_hat = mtk_estimate_channel(method, params)
    switch method
        case 'ls'
            H = params.H;
            Phi = params.Phi;
            N = params.N;
            rho = params.rho;

            Y = sqrt(rho)*H*conj(Phi') + N;
            R = mtk_util_quantize(Y, 1/sqrt(2));

            H_hat = (1/sqrt(rho))*(R*conj(Phi))/(conj(Phi')*conj(Phi));
        case 'blmmse'
            H = params.H;
            Phi = params.Phi;
            N = params.N;
            rho = params.rho;
            C_h_real = params.C_h_real;

            M = length(H(:,1));
            K = length(H(1,:));
            tau = length(Phi(:,1));

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(2/(K*rho + 1)) * eye(2*M*tau);
            Phi_tilde_real = A_real * Phi_line_real;

            C_y_r = Phi_line_real * C_h_real * Phi_line_real' + 1/2 * eye(2*M*tau);
            Sigma_y_real = diag(diag(C_y_r));
            invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
            C_r_real = real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));

            Y = sqrt(rho)*H*conj(Phi') + N;
            
            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(y_real, 1);
            
            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat = reshape(h_hat, [M, K]);
        case 'additive-quantizer-noise'
            H = params.H;
            Phi = params.Phi;
            N = params.N;
            rho = params.rho;
            C_h_real = params.C_h_real;

            M = length(H(:,1));
            K = length(H(1,:));
            tau = length(Phi(:,1));

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(2/(K*rho + 1)) * eye(2*M*tau);
            Phi_tilde_real = A_real * Phi_line_real;

            C_y_r = Phi_line_real * C_h_real * Phi_line_real' + 1/2 * eye(2*M*tau);
            % Sigma_y_real = diag(diag(C_y_r));
            % invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
            % C_r_real = real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));
            C_r_real = (1-2/pi)*eye(2*tau*M) + A_real*C_y_r*A_real';

            Y = sqrt(rho)*H*conj(Phi') + N;
            
            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(y_real, 1);
            
            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat = reshape(h_hat, [M, K]);
        otherwise
            error('Channel type not implemented');
    end
end

