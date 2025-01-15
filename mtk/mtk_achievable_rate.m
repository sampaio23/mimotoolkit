function I = mtk_achievable_rate(params)
    B = params.B;    
    H = params.H;
    H_hat = params.H_hat;
    K = params.K;
    rho = params.rho;
    
    
    Hr = [real(H) -imag(H); imag(H) real(H)];
    Hr_hat = [real(H_hat) -imag(H_hat); imag(H_hat) real(H_hat)];
    
    C_zr = 1/2 * rho * B * (Hr * Hr') * B' + 1/2 * (B*B');
    Kr = sqrt(inv(diag(diag(C_zr))));
    C_zq = 2/pi * real(asin(Kr * C_zr * Kr));
    C_zq_xr = sqrt(2/pi)*rho/2*Kr*B*Hr;
    G_r = C_zq\C_zq_xr;
    A_r_d = sqrt(2/pi)*Kr;      
    C_n_r_q_d = C_zq - A_r_d * C_zr * A_r_d';
    E_r = Hr - Hr_hat;
    
    I = 0;
    for k=1:2*K
        g_r_k = G_r(:,k);
        d_r_k = g_r_k' * A_r_d * B;
        h_hat_r_k = Hr_hat(:,k);
        
        num = rho*norm(d_r_k*h_hat_r_k)^2;
        den = norm(d_r_k)^2 + 2*g_r_k'*C_n_r_q_d*g_r_k;
        for i=1:2*K
            e_r_i = E_r(:,i);
            h_hat_r_i = Hr_hat(:,i);
            den = den + rho*norm(d_r_k * e_r_i)^2;
            if i ~= k
                den = den + rho*norm(d_r_k * h_hat_r_i)^2;
            end
        end

        I = I + 1/2*log2(1 + (num)/(den));
    end
end