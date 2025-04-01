
load dollarkurs.mat
X = USDSEK;
N = length(X);
tt=(1:N)';

%% 3a Linjär modell
A = [ones(N,1),tt];

mkdum = (A'*A) \ (A'*X);
mksmart = A \ X;

disp("Skillnader i koefficienter");
disp(mkdum-mksmart);

disp("Konditionstal: ")
disp(cond(A));

Y=tt*mkdum(2) + mkdum(1);
figure(1);
scatter(tt, X);
hold on;
plot(tt, Y);
errs = Y - X;
plot(tt, errs);
hold off;

%% 3b Normalekvationerna och konditionering

%% 3c Linjär + periodisk modell

% Er kod här...

L = 30;
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
function [Y] = F(d0, d1, d2, d3, L)
    load dollarkurs.mat
    X = USDSEK;
    N = length(X);
    tt=(1:N)';
    Y=zeros(N,1);
    for i = 1:N
        Y(i)=d0+d1*i+d2*sin(2*pi*i/L)+d3*cos(2*pi*i/L)-X(i);
    end
end

function [J] = jac(d0,d1,d2,d3,L)
    load dollarkurs.mat
    X = USDSEK;
    N = length(X);
    tt=(1:N)';
    Y=zeros(N,1);

    J = zeros(N, 5);
    fprintf("%f,%f,%f,%f,%f\n", d0,d1,d2,d3,L);

    for i = 1:N
        disp(i);
        J(i,1)=1;
        J(i,2)=i;
        disp(sin(2*pi*i/L));
        J(i,3)=sin(2*pi*i/L);
        J(i,4)=cos(2*pi*i/L);
        J(i,5)=-d2*cos(2*pi*i/L)*2*pi*i/L^2 + d3*2*pi*i/L^2*sin(2*pi*i/L);
    end
end

S=[mkdum;L];
dxnorm = inf;
tau = 1e-5;
while (dxnorm>tau && iter<1000)
    dx = -jac(S(1),S(2),S(3),S(4),S(5))\F(S(1),S(2),S(3),S(4),S(5)); %lös med minsta kvadrat i varje steg
    dxnorm = norm(dx);
    S = S + dx;
    iter = iter+1;
    %disp([iter dxnorm])
end

figure(3);
scatter(tt, X);
hold on;
Yperiodgood = S(1) + S(2)*tt + S(3)*sin(2*pi*tt/L) + S(4)*cos(2*pi*tt/L);
plot(tt,Yperiodgood);
Yperiodbad = mkdum(1) + mkdum(2)*tt + mkdum(3)*sin(2*pi*tt/L) + mkdum(4)*cos(2*pi*tt/L);
plot(tt,Yperiodbad);
Ylinregged=tt*mkdum(2) + mkdum(1);
plot(tt,Ylinregged);
hold off;

disp("koefficienter");
disp(S);

disp("medelkvadratfel");
SqErrPeriodGood = sum((Yperiodgood - X).^2) / length(Yperiodgood);
SqErrPeriodBad = sum((Yperiodbad - X).^2) / length(Yperiodbad);
SqErrLinreg = sum((Ylinregged - X).^2) / length(Ylinregged);

fprintf("för bra periodisk: %.16f\n", SqErrPeriodGood);
fprintf("för dålig periodisk: %.16f\n", SqErrPeriodBad);
fprintf("för linjär: %.16f\n", SqErrLinreg);
