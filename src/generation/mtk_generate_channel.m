function [H, results] = mtk_generate_channel(model, params, seed)
    K = params.K;
    M = params.M;
    channel_iterations = params.channel_iterations;
    switch model
        case 'rayleigh'
            rng(seed);
            H = 1/sqrt(2)*(randn(M, K, channel_iterations) + 1i*randn(M, K, channel_iterations));
            results.C_h = eye(M*K);
            results.C_h_real = 1/2*eye(2*M*K);
            results.C_H = params.M * eye(params.K);
        case 'rician'
            rng(seed);
            H = zeros(M, K, channel_iterations);
            K_factor = params.K_factor;
            for it=1:channel_iterations
                for k=1:K
                    mu = sqrt(K_factor(k)/(2*(K_factor(k) + 1)));
                    sigma = sqrt(1/(2*(K_factor(k) + 1)));
                    gk = sigma*(randn(M, 1) + 1i*randn(M, 1)) + mu;
                    H(:, k, it) = gk;
                end
            end
            results.C_h_real = 1/2*eye(2*M*K);
        case 'large-scale'
            rng(seed);
            H = 1/sqrt(2)*(randn(M, K, channel_iterations) + 1i*randn(M, K, channel_iterations));
            for i=1:channel_iterations
                H(:,:,i) = H(:,:,i)*sqrt(diag(params.beta));
            end
            c = zeros(M*K, 1);
            for i=1:K
                c((i-1)*M+1:i*M) = kron(params.beta(i), ones(M, 1));
            end
            C = diag(c);
            results.C_h_real = 1/2 * mtk_util_mat_real(C);
        case 'kron-markov'
            r_k = params.r_k;
            channel_index = params.channel_index;
            H = zeros(M, K, channel_index, channel_iterations);
            eta = params.eta;

            etaM = kron(diag(eta), eye(params.M));
            results.eta_real = [real(etaM) -imag(etaM); imag(etaM) real(etaM)];
            zetas = sqrt(1-eta.^2);
            zetaM = kron(diag(zetas), eye(params.M));
            results.zeta_real = [real(zetaM) -imag(zetaM); imag(zetaM) real(zetaM)];

            rng(seed);
            R = zeros(M*K);
            for k=1:K
                theta = 2*pi*rand();
                Rk = ones(M);
                for i=1:M
                    for j=(i+1):M
                        Rk(i,j) = (r_k*exp( 1i*theta))^(j-i);
                        Rk(j,i) = (r_k*exp(-1i*theta))^(j-i);
                    end
                end
                R(M*(k-1)+1:M*k, M*(k-1)+1:M*k) = Rk;
            end
            results.C_h_real = 1/2 * [real(R) -imag(R); imag(R) real(R)];

            for it=1:channel_iterations
                for k=1:K
                    gk = 1/sqrt(2)*(randn(M, 1) + 1i*randn(M, 1));
                    H(:, k, 1, it) = sqrtm(R(M*(k-1)+1:M*k, M*(k-1)+1:M*k)) * gk;
                end

                for idx=2:channel_index
                    H_inn = zeros(M, K);
                    for k=1:K
                        gk = 1/sqrt(2)*(randn(M, 1) + 1i*randn(M, 1));
                        H_inn(:, k) = sqrtm(R(M*(k-1)+1:M*k, M*(k-1)+1:M*k)) * gk;
                    end
                    h_inn = H_inn(:);
                    h = H(:, :, idx - 1, it);
                    h = h(:);
                    h = etaM * h + zetaM * h_inn;
                    H(:, :, idx, it) = reshape(h, [M, K]);
                end
            end
        otherwise
            error('Channel type not implemented');
    end
end
