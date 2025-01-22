function nmse = mtk_sim_ce_index(sim)
    method = sim.method;
    params = sim.params;
    params.quantizer = sim.quantizer;
    channel_iterations = params.channel_iterations;
    channel_index = params.channel_index;

    H = params.H;
    N = params.N;
    Phi = params.Phi;

    % TODO: move outside
    Beff = eye(2*params.M*params.tau);
    alpha = 0;

    params.rho = mtk_util_db_to_linear(params.SNR);
    params = mtk_ce_prepare(method, params);
    eta_real = params.eta_real;
    zeta_real = params.zeta_real;
    C_h_real = params.C_h_real;
    Cn_tilde_real = params.Cn_tilde_real;
    C_r_real = params.C_r_real;
    Phi_tilde_real = params.Phi_tilde_real;

    errors_final = zeros(channel_index, 1);
    Ms_final = zeros(channel_index, 1);

    errors_mean = zeros(channel_index, channel_iterations);
    Ms_mean = zeros(channel_index, channel_iterations);

    for it=1:channel_iterations
        h = H(:,:,1,it);
        h = h(:);
        h_real = [real(h); imag(h)];

        Y = sqrt(params.rho)*H(:,:,1,it)*conj(Phi') + N(:,:,1);

        y = Y(:);
        y_real = [real(y); imag(y)];

        r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

        h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
        Mk = C_h_real - C_h_real * (Phi_tilde_real' / C_r_real) * Phi_tilde_real * C_h_real;

        errors = zeros(channel_index, 1);
        Ms = zeros(channel_index, 1);

        errors(1) = norm(h_hat_real - h_real)^2 / (params.M * params.K);
        Ms(1) = trace(Mk) / (params.M * params.K);

        for i=2:channel_index
            h = H(:,:,i,it);
            h = h(:);
            h_real = [real(h); imag(h)];

            Y = sqrt(params.rho)*H(:,:,i,it)*conj(Phi') + N(:,:,i);

            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

            % Kalman
            h_hati_real = eta_real * h_hat_real;
            Mi = eta_real * Mk * eta_real' + zeta_real * C_h_real * zeta_real';
            Ki = Mi * Phi_tilde_real' / (Cn_tilde_real + Phi_tilde_real * Mi * Phi_tilde_real');
            h_hat_real = h_hati_real + Ki * (r - Phi_tilde_real * h_hati_real);
            Mk = (eye(2*params.M*params.K) - Ki * Phi_tilde_real) * Mi;

            errors(i) = norm(h_hat_real - h_real)^2 / (params.M * params.K);
            Ms(i) = trace(Mk) / (params.M * params.K);
        end

        errors_mean(:, it) = errors;
        Ms_mean(:, it) = Ms;

    end

    for i=1:channel_index
        errors_final(i) = mean(errors_mean(i,:));
        Ms_final(i) = mean(Ms_mean(i,:));
    end

    % FIXME: test
    nmse = Ms_final;
end

