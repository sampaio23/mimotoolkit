function [B, Beff] = mtk_util_comparator_network(type, params, seed)
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

            Brp = C(:,1:M);
            Bip = C(:,M+1:2*M);
            Bp = [Brp Bip];
            B = [eye(2*M); 1/sqrt(2)*Bp];
            Beffrp = kron(Brp, eye(tau));
            Beffip = kron(Bip, eye(tau));
            Beffp = [Beffrp Beffip];
            Beff = [eye(tau*2*M); 1/sqrt(2)*Beffp];
        case 'random'
            F = zeros(n, 2*M);
            i = 1;
            ii = 2;
            start_j = 2;
            for row=1:n
               F(row,i) = 1;
               F(row,ii) = -1;
               if (ii == 2*M)
                    i = i+1;
                    start_j = start_j+1;
                    ii = start_j;
               else
                   ii = ii+1;
               end
            end

            rng(seed);
            rows = randperm(n,alpha);

            rows = sort(rows);
            for i=1:length(rows)
                C(i,:) = F(rows(i),:);
            end

            Brp = C(:,1:M);
            Bip = C(:,M+1:2*M);
            Bp = [Brp Bip];
            B = [eye(2*M); 1/sqrt(2)*Bp];
            Beffrp = kron(Brp, eye(tau));
            Beffip = kron(Bip, eye(tau));
            Beffp = [Beffrp Beffip];
            Beff = [eye(tau*2*M); 1/sqrt(2)*Beffp];
        case 'dynamic'
            F = zeros(n, 2*M);
            i = 1;
            ii = 2;
            start_j = 2;
            for row=1:n
               F(row,i) = 1;
               F(row,ii) = -1;
               if (ii == 2*M)
                    i = i+1;
                    start_j = start_j+1;
                    ii = start_j;
               else
                   ii = ii+1;
               end
            end

            Beff = [eye(tau*2*M); zeros(tau*alpha, tau*2*M)];

            rng(seed);

            for t=1:tau
                rows = randperm(n,alpha);
    
                rows = sort(rows);
                for i=1:length(rows)
                    C(i,:) = F(rows(i),:);
                end

                current_tau = zeros(1,tau);
                current_tau(t) = 1;

                for j=1:length(rows)
                    Beff(tau*2*M + (j-1)*tau + t,:) = 1/sqrt(2)*kron(C(j,:), current_tau);
                end
            end

            B = []; % Dynamic only used for channel estimation
        otherwise
            error('Comparator network type not implemented');
    end

    
end
