function phi = mtk_generate_pilot(type, params)
    switch type
        case 'dft'
            K = params.K;
            tau = params.tau;

            tmp = dftmtx(tau);
            phi = tmp(:, 1:K);
        otherwise
            error('Pilot type not implemented');
    end
end

