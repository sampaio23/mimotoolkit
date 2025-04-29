function phi = mtk_generate_pilot(type, params)
    K = params.K;
    tau = params.tau;    
    switch type
        case 'dft'
            tmp = dftmtx(tau);
            phi = tmp(:, 1:K);
        case 'gaussian'
            rng(0);
            phi = 1/sqrt(2)*(rand(tau, K) + 1i*rand(tau, K));
        otherwise
            error('Pilot type not implemented');
    end
end

