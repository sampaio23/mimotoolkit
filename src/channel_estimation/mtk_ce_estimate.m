function [H_hat, results] = mtk_ce_estimate(method, params)
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
        case 'kfb'
            C_h_real = params.C_h_real;
            H = params.H;
            K = params.K;
            M = params.M;
            N = params.N;
            Phi = params.Phi;

            eta_real = params.eta_real;
            zeta_real = params.zeta_real;
            Phi_tilde_real = params.Phi_tilde_real;
            Cn_tilde_real = params.Cn_tilde_real;
            C_r_real = params.C_r_real;

            % TODO: remove and unify
            Beff = eye(2*params.M*params.tau);
            alpha = 0;

            H_hat = zeros(M, K, params.channel_index);
            results.mmse = zeros(2*M*K, 2*M*K, params.channel_index);

            Y = sqrt(params.rho)*H(:,:,1)*conj(Phi') + N(:,:,1);

            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

            Mk = C_h_real - C_h_real * (Phi_tilde_real' / C_r_real) * Phi_tilde_real * C_h_real;
            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat(:,:,1) = reshape(h_hat, [M, K]);
            results.mmse(:,:,1) = Mk;

            for i=2:params.channel_index
                Y = sqrt(params.rho)*H(:,:,i)*conj(Phi') + N(:,:,i);

                y = Y(:);
                y_real = [real(y); imag(y)];

                r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

                % Kalman
                h_hati_real = eta_real * h_hat_real;
                Mi = eta_real * Mk * eta_real' + zeta_real * C_h_real * zeta_real';
                Ki = Mi * Phi_tilde_real' / (Cn_tilde_real + Phi_tilde_real * Mi * Phi_tilde_real');
                h_hat_real = h_hati_real + Ki * (r - Phi_tilde_real * h_hati_real);
                Mk = (eye(2*params.M*params.K) - Ki * Phi_tilde_real) * Mi;

                h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

                H_hat(:,:,i) = reshape(h_hat, [M, K]);
                results.mmse(:,:,i) = Mk;
            end
        case 'kfb-cn'
            C_h_real = params.C_h_real;
            H = params.H;
            K = params.K;
            M = params.M;
            N = params.N;
            Phi = params.Phi;

            eta_real = params.eta_real;
            zeta_real = params.zeta_real;
            Phi_tilde_real = params.Phi_tilde_real;
            Cn_tilde_real = params.Cn_tilde_real;
            C_r_real = params.C_r_real;

            Beff = params.Beff;
            alpha = params.alpha;

            H_hat = zeros(M, K, params.channel_index);
            results.mmse = zeros(2*M*K, 2*M*K, params.channel_index);

            Y = sqrt(params.rho)*H(:,:,1)*conj(Phi') + N(:,:,1);

            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

            Mk = C_h_real - C_h_real * (Phi_tilde_real' / C_r_real) * Phi_tilde_real * C_h_real;
            h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
            h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

            H_hat(:,:,1) = reshape(h_hat, [M, K]);
            results.mmse(:,:,1) = Mk;

            for i=2:params.channel_index
                Y = sqrt(params.rho)*H(:,:,i)*conj(Phi') + N(:,:,i);

                y = Y(:);
                y_real = [real(y); imag(y)];

                r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

                % Kalman
                h_hati_real = eta_real * h_hat_real;
                Mi = eta_real * Mk * eta_real' + zeta_real * C_h_real * zeta_real';
                Ki = Mi * Phi_tilde_real' / (Cn_tilde_real + Phi_tilde_real * Mi * Phi_tilde_real');
                h_hat_real = h_hati_real + Ki * (r - Phi_tilde_real * h_hati_real);
                Mk = (eye(2*params.M*params.K) - Ki * Phi_tilde_real) * Mi;

                h_hat = h_hat_real(1 : M*K) + 1i * h_hat_real(M*K+1 : 2*M*K);

                H_hat(:,:,i) = reshape(h_hat, [M, K]);
                results.mmse(:,:,i) = Mk;
            end
        otherwise
            error('Channel type not implemented');
    end
end

