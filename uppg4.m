

%% 4a Konvergensstudie
% Er kod här...
NN=100;
Ns = round(logspace(1, 3, NN));
Is = ones(NN,3);
for i=1:NN
    Is(i,1)=trapets(0.11,Ns(i));
    Is(i,2)=trapets(0.32, Ns(i));
    Is(i,3)=trapets(1.14,Ns(i));
end
erf1=erf(0.11);
erf2=erf(0.32);
erf3=erf(1.14);

figure(1)
loglog(Ns, abs(Is(:,1) - erf1), "--", "Color",[0.3010 0.7450 0.9330]);
hold on;
loglog(Ns, abs(Is(:,2) - erf2), "--", "Color",[0.8500 0.3250 0.0980]);
loglog(Ns, abs(Is(:,3) - erf3), "--", "Color",[0.9290 0.6940 0.1250]);

C= 1/(3*sqrt(pi));
f=@(x, N) C * x^3 ./ (N .^ 2);
loglog(Ns, f(0.11, Ns), "Color",[0.3010 0.7450 0.9330]);
loglog(Ns, f(0.32, Ns), "Color",[0.8500 0.3250 0.0980]);
loglog(Ns, f(1.14, Ns), "Color",[0.9290 0.6940 0.1250]);


hold off

Xs = linspace(0,6,NN);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I=trapets(x, N)
%  Indata:
%
%  N  - antal delintervall (skalär)
%
%  Utdata:
%
%  I - integralen (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=linspace(0,x,N+1);
    g = @(x) 2/sqrt(pi)*exp(-(x.^2));
    Y = g(X);
    y1 = Y(1)/2;
    yn = Y(end)/2;

    I = x/N*(y1+yn + sum(Y(2:end-1)));
end








