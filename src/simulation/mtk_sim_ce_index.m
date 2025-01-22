function [nmse, results] = mtk_sim_ce_index(sim)
    method = sim.method;
    params = sim.params;
    params.quantizer = sim.quantizer;
    channel_index = params.channel_index;

    H = params.H;
    N = params.N;

    params.rho = mtk_util_db_to_linear(params.SNR);

    params = mtk_ce_prepare(method, params);

    nmse = zeros(channel_index, 1);
    results.mmse = zeros(channel_index, 1);
    for it=1:params.channel_iterations
        params.H = H(:,:,:,it);
        params.N = N(:,:,(it-1)*channel_index+1:it*channel_index);

        [H_hat, results_ce] = mtk_ce_estimate(method, params);

        for i=1:channel_index
            E = params.H(:,:,i) - H_hat(:,:,i);
            nmse(i) = nmse(i) + norm(E, 'fro')^2;
            results.mmse(i) = results.mmse(i) + trace(results_ce.mmse(:,:,i));
        end
    end
    nmse = nmse / (params.channel_iterations * params.K * params.M);
    results.mmse = results.mmse / (params.channel_iterations * params.K * params.M);
end

