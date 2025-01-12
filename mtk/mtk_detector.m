function G = mtk_detector(detector, params)
    H_hat = params.H_hat;    
    switch detector
        case 'zf'        
            G = pinv(H_hat);
        case {'mf', 'mrc'}
            G = inv(diag(diag(H_hat' * H_hat))) * H_hat';
        case 'mmse'
            G = (1/params.rho*eye(params.K) + H_hat'*H_hat)\H_hat';
        case 'bmrc'
            S_r = (H_hat * H_hat') + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = inv(diag(diag(A' * A))) * A';
        case 'blmmse'
            S_r = (H_hat * H_hat') + 1/params.rho*eye(params.M);
            S_y = sqrt(inv(diag(diag(S_r))));
            A = sqrt(2/pi) * S_y * H_hat;
            G = A'/(A*A' + 2/pi*(asin(S_y * S_r * S_y) - S_y * S_r * S_y + 1/params.rho * inv(diag(diag(S_r)))));
        otherwise
            error('Detector not implemented');
    end
end

