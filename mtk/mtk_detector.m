function G = mtk_detector(detector, params)
    H_hat = params.H_hat;    
    switch detector
        case 'zf'        
            G = pinv(H_hat);
        case {'mf', 'mrc'}
            G = H_hat';
        case 'mmse'
            G = (1/params.rho*eye(params.K) + H_hat'*H_hat)\H_hat';
        case 'blmmse'
            C_z = params.rho * H_hat * H_hat' + eye(params.M);
            S_y = sqrt(inv(diag(diag(C_z))));
            G = sqrt(pi/2)*sqrt(params.rho)*H_hat'*S_y/(asin(S_y * real(C_z) * S_y) + 1i*asin(S_y * imag(C_z) * S_y));
        otherwise
            error('Detector not implemented');
    end
end

