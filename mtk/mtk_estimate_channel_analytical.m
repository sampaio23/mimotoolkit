function nmse = mtk_estimate_channel_analytical(method, params)
    switch method
        case {'blmmse', 'blmmse-cn'}
            C_h_real = params.C_h_real;
            C_r_real = params.C_r_real;
            K = params.K;
            M = params.M;
            Phi_tilde_real = params.Phi_tilde_real;
            nmse = 1/(M*K)*trace(C_h_real - C_h_real * (Phi_tilde_real' / C_r_real) * Phi_tilde_real * C_h_real);
        otherwise
            error('Channel type not implemented');
    end
end

