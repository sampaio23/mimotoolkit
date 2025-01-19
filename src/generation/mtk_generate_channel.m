function H = mtk_generate_channel(model, params, seed)
    K = params.K;
    M = params.M;
    channel_iterations = params.channel_iterations;
    switch model
        case 'rayleigh'
            rng(seed);
            H = 1/sqrt(2)*(randn(M, K, channel_iterations) + 1i*randn(M, K, channel_iterations));
        case 'kron-markov'
            r_k = params.r_k;
            channel_index = params.channel_index;
            H = zeros(M, K, channel_index, channel_iterations);
            eta = params.eta;

            etaM = kron(diag(eta), eye(M));
            zetas = sqrt(1-eta.^2);
            zetaM = kron(diag(zetas), eye(M));

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
