function nmse = mtk_sim_ce_index(sim)
    method = sim.method;
    params = sim.params;
    params.quantizer = sim.quantizer;
    channel_iterations = params.channel_iterations;
    channel_index = params.channel_index;

    H = params.H;
    N = params.N;

    % TODO: move to generation, return result struct
    eta = params.eta;
    etaM = kron(diag(eta), eye(params.M));
    etar = [real(etaM) -imag(etaM); imag(etaM) real(etaM)];
    zetas = sqrt(1-eta.^2);
    zetaM = kron(diag(zetas), eye(params.M));
    zetar = [real(zetaM) -imag(zetaM); imag(zetaM) real(zetaM)];


    K = params.K;
    M = params.M;
    tau = params.tau;


    Phi = params.Phi;

    Beff = eye(2*params.M*params.tau);
    alpha = 0;

    rho = mtk_util_db_to_linear(params.SNR);

    %Phi_tilde_real = params.Phi_tilde_real;
    %C_r_real = params.C_r_real;
    C_h_real = params.C_h_real; % C_h_real = Rr/2;

    Phi_line = kron(Phi, sqrt(rho)*eye(M));
    Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
    A_real = sqrt(2/pi) * sqrt(1/(K*rho + 1)) * eye((2*M + alpha)*tau);
    Phi_tilde_real = A_real * Beff * Phi_line_real;

    C_y_r = Beff * (Phi_line_real * C_h_real * Phi_line_real' + 1/2*eye(2*M*tau)) * Beff';
    Sigma_y_real = diag(diag(C_y_r));
    invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
    C_r_real = 1/2 * real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));
    Cn_tilde_real = C_r_real - Phi_tilde_real * C_h_real * Phi_tilde_real';
    % The 1/2 factor comes from 1/sqrt(2) quantization factor (?)

    errors_final = zeros(channel_index, 1);
    Ms_final = zeros(channel_index, 1);

    errors_mean = zeros(channel_index, channel_iterations);
    Ms_mean = zeros(channel_index, channel_iterations);


    for it=1:channel_iterations
        h = H(:,:,1,it);
        h = h(:);
        h_real = [real(h); imag(h)];

        Y = sqrt(rho)*H(:,:,1,it)*conj(Phi') + N(:,:,1);

        y = Y(:);
        y_real = [real(y); imag(y)];

        r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

        h_hat_real = C_h_real * (Phi_tilde_real' / C_r_real) * r;
        Mk = C_h_real - C_h_real * (Phi_tilde_real' / C_r_real) * Phi_tilde_real * C_h_real;

        errors = zeros(channel_index, 1);
        Ms = zeros(channel_index, 1);

        errors(1) = norm(h_hat_real - h_real)^2 / (M*K);
        Ms(1) = trace(Mk) / (M*K);

        for i=2:channel_index
            h = H(:,:,i,it);
            h = h(:);
            h_real = [real(h); imag(h)];

            Y = sqrt(rho)*H(:,:,i,it)*conj(Phi') + N(:,:,i);

            y = Y(:);
            y_real = [real(y); imag(y)];

            r = mtk_util_quantize(params.quantizer, Beff * y_real, 1/sqrt(2));

            % Kalman
            h_hati_real = etar * h_hat_real;
            Mi = etar * Mk * etar' + zetar * C_h_real * zetar';
            Ki = Mi * Phi_tilde_real' / (Cn_tilde_real + Phi_tilde_real * Mi * Phi_tilde_real');
            h_hat_real = h_hati_real + Ki * (r - Phi_tilde_real * h_hati_real);
            Mk = (eye(2*M*K) - Ki * Phi_tilde_real) * Mi;

            errors(i) = norm(h_hat_real - h_real)^2 / (M*K);
            Ms(i) = trace(Mk) / (M*K);
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

