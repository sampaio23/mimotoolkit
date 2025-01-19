function SER = mtk_sim_detect_snr(sim)
    method = sim.method;
    modulation = sim.modulation;
    detector = sim.detector;
    quantizer = sim.quantizer;
    SNRs = sim.SNRs;
    params = sim.params;

    H = params.H;
    n = params.n;

    SER = zeros(size(SNRs));
    snr_index = 1;
    for SNR = SNRs
        current_error_bmmse = 0;
    
        params.rho = mtk_util_db_to_linear(SNR);
        for channel_index=1:params.channel_iterations
            params.H = H(:,:,channel_index);
            params.H_hat = mtk_ce_estimate(method, params);
    
            G = mtk_det_detector(detector, params);
            for noise_index=1:params.noise_iterations
                params.n = repmat(n(:,:,noise_index), 1, params.symbol_iterations);

                r = sqrt(params.rho) * params.H * params.x + params.n;
                y = mtk_util_quantize(quantizer, r, (params.K*params.rho + 1)/(pi*params.rho));
    
                x_tilde = G * y;
                x_hat = mtk_util_demodulate(modulation, x_tilde);
                current_error_bmmse = current_error_bmmse + sum(x_hat ~= params.x, 'all');
            end
        end
        SER(snr_index) = current_error_bmmse / (params.K * params.channel_iterations * params.symbol_iterations * params.noise_iterations);
        snr_index = snr_index + 1;
    end
end

