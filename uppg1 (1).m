
%%
function [mu, iter] = potens(A,tau)
    [m, m] = size(A);
    A = sparse(A);
    init = ones(m);
    done = False;
    iter = 0;
    while not done
        init2 = A * init;
        mu = init * init2;
        init2 = init2 ./ norm(init2);
        iter = iter + 1;
        if norm(init2 - init) < tau
            done = True; % replace with some better tolerance check
        end
    end
end

function [mu, iter] = inverspotens(A, tau)
    % formulation should be equiv to regular power method with inverse matrix
    A = sparse(A);
    A_inv = inv(A);
    [mu_inv, iter] = potens(A_inv, tau);
    mu = 1 / mu_inv;
end

%% 1d -- Beräkning av andra egenvärden

function [mu, iter] = inverspotens_skift(A, skift, tau)
    % reuses inverse power method again
    I = speye(size(A));
    [mu_shifted, iter] = inverspotens(A - skift * I, tau);
    mu = skift + 1 / mu_shifted;
end

