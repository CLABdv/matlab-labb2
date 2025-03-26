%% 1b -- Visualisering

% found on matlab documentation
[V,D] = eig(A);
[d, ind] = sort(diag(D)); %sorted eigenvalues and permutationvector
Vsorted = V(:,ind); %vectors sorted by permutationvector


load("eiffel1.mat");
for i=[1,2,3,4,100]
    figure(i)
    e = Vsorted(:,i);
    trussplot(xnod+e(1:2:end),ynod+e(2:2:end),bars)
end

%% 1c -- Beräkning av största och minsta egenvärdena

% Er kod här...


%% 1d -- Beräkning av andra egenvärden

% Er kod här...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = potens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (skälär)
%
%  Utdata:
%
%  mu - största egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Er kod här...

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = inverspotens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (skälär)
%
%  Utdata:
%
%  mu - minsta egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Er kod här...

end
