

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

%%
% Kod för approximation av noggrannhetsordningen
n = 15;
x = 0.11;
a = linspace(1, n, n);
a = 2 .^ a; % dubbla finheten i integrationen för varje steg
res = zeros(2, n-2);
for i = 1:n-2
    res(1, i) = a(i);
    res(2,i) = (trapets(x, a(i)) - trapets(x, a(i+1))) / (trapets(x, a(i+1)) - trapets(x, a(i + 2)));
end
T = array2table(res);
disp(T)
%%
% Kod för att se noggrannheten variera i x med fixat N
points_eval = linspace(0, 6, 100); % punkterna där erf evalueras
errs = zeros(3, 100);
Ns = [50, 120, 400];
for i = 1:100
    for j = 1:3
        errs(j,i) = abs(trapets(points_eval(i), Ns(j)) - erf(points_eval(i)));
    end
end

figure(2)
semilogy(points_eval, errs(1,:))
hold on
semilogy(points_eval, errs(2,:))
semilogy(points_eval, errs(3,:))


%%




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








