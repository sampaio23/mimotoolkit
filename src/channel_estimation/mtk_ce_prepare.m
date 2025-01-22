function params = mtk_ce_prepare(method, params)
    switch method
        case 'ls'
        case 'blmmse'
            C_h_real = params.C_h_real;
            K = params.K;
            M = params.M;
            Phi = params.Phi;
            rho = params.rho;
            tau = params.tau;

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(2/(K*rho + 1)) * eye(2*M*tau);
            params.Phi_tilde_real = A_real * Phi_line_real;

            C_y_r = Phi_line_real * C_h_real * Phi_line_real' + 1/2 * eye(2*M*tau);
            Sigma_y_real = diag(diag(C_y_r));
            invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
            params.C_r_real = real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));
        case 'aqnm'
            C_h_real = params.C_h_real;
            K = params.K;
            M = params.M;
            Phi = params.Phi;
            rho = params.rho;
            tau = params.tau;

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(2/(K*rho + 1)) * eye(2*M*tau);
            params.Phi_tilde_real = A_real * Phi_line_real;

            C_y_r = Phi_line_real * C_h_real * Phi_line_real' + 1/2 * eye(2*M*tau);
            params.C_r_real = (1-2/pi)*eye(2*tau*M) + A_real*C_y_r*A_real';
        case 'blmmse-cn'
            alpha = params.alpha;
            Beff = params.Beff;
            C_h_real = params.C_h_real;
            K = params.K;
            M = params.M;
            Phi = params.Phi;
            rho = params.rho;
            tau = params.tau;

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(2/(K*rho + 1)) * eye((2*M + alpha) * tau);
            params.Phi_tilde_real = A_real * Beff * Phi_line_real;

            C_y_r = Beff * (Phi_line_real * C_h_real * Phi_line_real' + 1/2 * eye(2*M*tau)) * Beff';
            Sigma_y_real = diag(diag(C_y_r));
            invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
            params.C_r_real = real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));
        case 'kfb'
            C_h_real = params.C_h_real;
            K = params.K;
            M = params.M;
            Phi = params.Phi;
            rho = params.rho;
            tau = params.tau;

            % TODO: move outside
            Beff = eye(2*params.M*params.tau);
            alpha = 0;

            Phi_line = kron(Phi, sqrt(rho)*eye(M));
            Phi_line_real = [ real(Phi_line) -imag(Phi_line); imag(Phi_line) real(Phi_line)];
            A_real = sqrt(2/pi) * sqrt(1/(K*rho + 1)) * eye((2*M + alpha)*tau);
            Phi_tilde_real = A_real * Beff * Phi_line_real;

            C_y_r = Beff * (Phi_line_real * C_h_real * Phi_line_real' + 1/2*eye(2*M*tau)) * Beff';
            Sigma_y_real = diag(diag(C_y_r));
            invsqr_Sigma_y_real = sqrt(inv(Sigma_y_real));
            C_r_real = 1/2 * real(2/pi * asin(invsqr_Sigma_y_real * real(C_y_r) * invsqr_Sigma_y_real));

            params.Cn_tilde_real = C_r_real - Phi_tilde_real * C_h_real * Phi_tilde_real';
            params.Phi_tilde_real = Phi_tilde_real;
            params.C_r_real = C_r_real;
        otherwise
            error('Channel type not implemented');
    end
end

