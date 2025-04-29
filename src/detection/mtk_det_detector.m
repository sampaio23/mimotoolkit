function G = mtk_det_detector(detector, params)
    H_hat = params.H_hat;    
    switch detector
        case 'zf'        
            G = pinv(H_hat);
        case {'mf', 'mrc'}
            G = (diag(diag(H_hat' * H_hat))) \ H_hat';
        case 'mmse'
            G = (1/params.rho*eye(params.K) + H_hat'*H_hat)\H_hat';
        case 'bzf'
            S_r = (H_hat * H_hat') + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = pinv(A);
        case 'bmrc'
            S_r = (H_hat * H_hat') + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = (diag(diag(A' * A))) \ A';
        case 'bmmse'
            S_r = (H_hat * H_hat') + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = A'/(2/pi*(asin(S_y * real(S_r) * S_y) + 1i*asin(S_y * imag(S_r) * S_y)));
        case 'robust-bmmse'
            S_r = (H_hat * H_hat') + params.Gamma + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = A'/(2/pi*(asin(S_y * real(S_r) * S_y) + 1i*asin(S_y * imag(S_r) * S_y)));
        case 'bmmse-cn'
            B = params.B;

            H_hat_real = mtk_util_mat_real(H_hat);
            S_r_real = 1/2 * (B * (H_hat_real * H_hat_real') * B') + 1/2 * 1/params.rho * B*B';
            S_y_real = sqrt(inv(diag(diag(S_r_real))));
            A = sqrt(2/pi) * 1/2 * S_y_real * B * H_hat_real;
            G = A'/(2/pi*(asin(S_y_real * real(S_r_real) * S_y_real)));
        case 'robust-bmmse-cn'
            B = params.B;

            H_hat_real = mtk_util_mat_real(H_hat);
            S_r = 1/2 * B * (H_hat_real * H_hat_real') * B' + B * mtk_util_mat_real(params.Gamma) * B' + 1/2 * 1/params.rho * (B*B');
            S_y = sqrt(inv(diag(diag(S_r))));
            A_real = sqrt(2/pi) * 1/2 * S_y * B * H_hat_real;
            G = A_real'/(2/pi*(asin(S_y * real(S_r) * S_y)));
        otherwise
            error('Detector not implemented');
    end
end

