
load dollarkurs.mat
X = USDSEK;
N = length(X);
tt=(1:N)';

%% 3a Linjär modell
A = [ones(N,1),tt];

mkdum = (A'*A) \ (A'*X);
mksmart = A \ X;



Ylinregged=tt*mkdum(2) + mkdum(1);
figure(1);
scatter(tt, X);
hold on;
plot(tt, Ylinregged);
errs = Ylinregged - X;
plot(tt, errs);
hold off;

%% 3b Normalekvationerna och konditionering

disp("Skillnader i koefficienter");
disp(mkdum-mksmart);

disp("Konditionstal: ")
disp(cond(A));

%% 3c Linjär + periodisk modell

% Er kod här...

L = 350;
A=[ones(N,1), tt, sin(2*pi*tt/L), cos(2*pi*tt/L)];

mkdum = (A'*A) \ (A'*X);
mksmart = A \ X;

disp("Skillnader i koefficienter");
disp(mkdum-mksmart);

disp("Konditionstal: ")
disp(cond(A));

Y = mkdum(1) + mkdum(2)*tt + mkdum(3)*sin(2*pi*tt/L) + mkdum(4)*cos(2*pi*tt/L);
figure(2);
scatter(tt, X);
hold on;
plot(tt,Y);
errs = Y - X;
%plot(tt,errs);
hold off;


%% 3d Icke-linjär modell
function [Y] = F(X)
    d0 = X(1);
    d1= X(2);
    d2 = X(3);
    d3 = X(4);
    L = X(5);
    load dollarkurs.mat
    X = USDSEK;
    N = length(X);
    Y=zeros(N,1);
    % TODO: Vectorise
    for i = 1:N
        Y(i)=d0+d1*i+d2*sin(2*pi*i/L)+d3*cos(2*pi*i/L)-X(i);
    end
end

function [J] = jac(X)
    d0 = X(1);
    d1= X(2);
    d2 = X(3);
    d3 = X(4);
    L = X(5);
    load dollarkurs.mat
    X = USDSEK;
    N = length(X);

    J = zeros(N, 5);
    fprintf("%f,%f,%f,%f,%f\n", d0,d1,d2,d3,L);

    for i = 1:N
        J(i,1)=1;
        J(i,2)=i;
        J(i,3)=sin(2*pi*i/L);
        J(i,4)=cos(2*pi*i/L);
        J(i,5)=2*pi*i/L^2*(-d2*cos(2*pi*i/L) + d3*sin(2*pi*i/L));
    end
end

S=[mkdum;L];
tau = 1e-13;
iter = 0;
while (true)
    dx = -jac(S) \ F(S); %lös med minsta kvadrat i varje steg
    if all (abs(dx) - tau < 0)
        break
    end
    S = S + dx;
    iter = iter+1;
   disp([iter dxnorm])
end

figure(3);
scatter(tt, X);
hold on;
v = [ones(N,1), tt, sin(2*pi*tt/L), cos(2*pi*tt/L)];

Yperiodgood = S(1) + S(2)*tt + S(3)*sin(2*pi*tt/S(5)) + S(4)*cos(2*pi*tt/S(5));
plot(tt,Yperiodgood);

Yperiodbad = sum(mkdum' .* v, 2);
plot(tt,Yperiodbad);

plot(tt,Ylinregged);
hold off;

legend({'y = USDSEK', 'y = Yperiodgood','y = Yperiodbad','Ylinregged'},'Location','northeast')

disp("koefficienter");
disp(S);

disp("medelkvadratfel");
SqErrPeriodGood = sum((Yperiodgood - X).^2) / length(Yperiodgood);
SqErrPeriodBad = sum((Yperiodbad - X).^2) / length(Yperiodbad);
SqErrLinreg = sum((Ylinregged - X).^2) / length(Ylinregged);

fprintf("för bra periodisk: %.16f\n", SqErrPeriodGood);
fprintf("för dålig periodisk: %.16f\n", SqErrPeriodBad);
fprintf("för linjär: %.16f\n", SqErrLinreg);
