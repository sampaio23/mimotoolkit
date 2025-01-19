function [B, Beff] = mtk_util_comparator_network(type, params)
    alpha = params.alpha;    
    M = params.M;
    tau = params.tau;

    n = nchoosek(2*M, 2); % Maximum number of comparators
    assert(alpha <= n, "You cannot choose more than the number of available comparators");

    C = zeros(alpha, 2*M);

    switch type
        case 'independent'
	        i = 1;
	        j = 2;
            for row=1:alpha
	           C(row,i) = 1;
	           C(row,j) = -1;
	           i = i+2;
	           j = j+2;
            end
        case 'random'
            F = zeros(n, 2*M);
            i = 1;
            ii = 2;
            start_j = 2;
            for row=1:n
               F(row,i) = 1;
               F(row,ii) = -1;
               if (ii == 2*Nr)
                    i = i+1;
                    start_j = start_j+1;
                    ii = start_j;
               else
                   ii = ii+1;
               end
            end

            current_seed = rng(9997*iterations);
            rows = randperm(n,alpha);
            rng(current_seed);

            rows = sort(rows);
            for i=1:length(rows)
                C(i,:) = F(rows(i),:);
            end
        otherwise
            error('Comparator network type not implemented');
    end

    Brp = C(:,1:M);
    Bip = C(:,M+1:2*M);
    Bp = [Brp Bip];
    B = [eye(2*M); 1/sqrt(2)*Bp];
    Beffrp = kron(Brp, eye(tau));
    Beffip = kron(Bip, eye(tau));
    Beffp = [Beffrp Beffip];
    Beff = [eye(tau*2*M); 1/sqrt(2)*Beffp];
end
