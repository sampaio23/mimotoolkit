function nmse = mtk_sim_ce_snr(sim)
    method = sim.method;
    SNRs = sim.SNRs;
    params = sim.params;
    params.quantizer = sim.quantizer;

    H = params.H;
    N = params.N;

    nmse = zeros(size(SNRs));
    snr_index = 1;
    for SNR = SNRs
        current_error_nmse = 0;
    
        params.rho = mtk_util_db_to_linear(SNR);
        params = mtk_ce_prepare(method, params);
        for channel_index=1:params.channel_iterations
            params.H = H(:,:,channel_index);
            params.N = N(:,:,1);
    
            H_hat = mtk_ce_estimate(method, params);
            E = H(:,:,channel_index) - H_hat;
            current_error_nmse = current_error_nmse + norm(E, 'fro')^2;
        end
        nmse(snr_index) = current_error_nmse/(params.channel_iterations * params.M * params.K);
        snr_index = snr_index + 1;
    end
end

