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

%% 1c -- Ber�kning av st�rsta och minsta egenv�rdena

% Er kod h�r...


%% 1d -- Ber�kning av andra egenv�rden

% Er kod h�r...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = potens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (sk�l�r)
%
%  Utdata:
%
%  mu - st�rsta egenv�rdet till A (skal�r)
%  iter - antal iterationer som anv�nts (skal�r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Er kod h�r...

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = inverspotens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (sk�l�r)
%
%  Utdata:
%
%  mu - minsta egenv�rdet till A (skal�r)
%  iter - antal iterationer som anv�nts (skal�r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Er kod h�r...

end
