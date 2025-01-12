function G = mtk_detector(detector, params)
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
            G = A'/(2/pi*(asin(S_y * S_r * S_y)));
        otherwise
            error('Detector not implemented');
    end
end

