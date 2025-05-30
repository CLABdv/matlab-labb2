%% 1b -- Visualisering

% found on matlab documentation
load("eiffel1.mat");
[V,D] = eig(A);
[d, ind] = sort(diag(D)); %sorted eigenvalues and permutationvector
Vsorted = V(:,ind); %vectors sorted by permutationvector

for i=[1,2,3,4,100]
    figure(i)
    e = Vsorted(:,i);
    trussplot(xnod+e(1:2:end),ynod+e(2:2:end),bars)
end

trussanim(xnod,ynod,bars,Vsorted(:,2));

%% 1c

function [mu, iter] = potens(A,tau)
    [m,~] = size(A);
    A = sparse(A);
    init = ones(m,1);
    done = false;
    iter = 0;
    mu2=NaN;
    while ~done
        init2 = A * init;
        mu=mu2;
        mu2 = init' * init2;
        init2 = init2 ./ norm(init2);
        iter = iter + 1;
        if abs(mu2-mu) < tau
            done = true; % replace with some better tolerance check
        end
        init=init2;
    end
end

function [mu, iter] = inverspotens(A, tau)
    % formulation should be equiv to regular power method with inverse matrix
    A = sparse(A);
    A_inv = inv(A);
    [mu_inv, iter] = potens(A_inv, tau);
    mu=1/mu_inv;
end

T=zeros(4,6);
for i=1:4
    file =sprintf("eiffel%d.mat",i);
    load(file);
    [~,D] = eig(A);
    d = sort(diag(D)); %sorted eigenvalues
    T(i,3) = d(2)/d(1);
    [T(i,1), T(i,2)] = potens(A,1e-10);
    [T(i,4), T(i,5)] = inverspotens(A,1e-10);
    T(i,6) = d(end)/d(end-1);
end

tab=array2table(T,'VariableNames',{'St�rsta egenv�rdet' '# iter' 'lambda_2/lambda_1' 'Minsta egenv�rdet', '#iter', 'lambda_n/lambda_{n-1}'},'RowNames',{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
disp(tab);
%% 1d -- Ber�kning av andra egenv�rden

function [mu, iter] = inverspotens_skift(A, skift, tau)
    % reuses inverse power method again
    I = speye(size(A));
    [mu_shifted, iter] = inverspotens(A - skift * I, tau);
    mu = skift + 1 / mu_shifted;
end
load("eiffel1.mat");
egen8 = inverspotens_skift(A,8,1e-5);
egen55 = inverspotens_skift(A,55,1e-5);
egen67 = inverspotens_skift(A,67,1e-5);

disp(egen8);
disp(egen55);
disp(egen67);