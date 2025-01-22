function nmse = mtk_sim_ce_index(sim)
    method = sim.method;
    params = sim.params;
    params.quantizer = sim.quantizer;
    channel_index = params.channel_index;

    H = params.H;
    N = params.N;

    params.rho = mtk_util_db_to_linear(params.SNR);

    params = mtk_ce_prepare(method, params);

    nmse = zeros(channel_index);
    for it=1:params.channel_iterations
        params.H = H(:,:,:,it);
        params.N = N;

        [H_hat, results] = mtk_ce_estimate(method, params);

        for i=1:channel_index
            E = params.H(:,:,i) - H_hat(:,:,i);
            % nmse(i) = nmse(i) + norm(E, 'fro')^2;
            nmse(i) = nmse(i) + trace(results.mmse(:,:,i));
        end
    end
    % nmse = nmse / (params.channel_iterations * params.K * params.M);
    nmse = nmse / (params.channel_iterations * params.K * params.M);
end

