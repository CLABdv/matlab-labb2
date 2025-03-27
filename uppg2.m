%% Interpolation

h = 0.25; % Får ej ändras i koden nedan

[t,x,y,vx,vy] = kastbana(h);
figure(1);
ps = [x,y];
scatter(x,y);

%% Linjär interpolation



x1 = ps(1,1);
y1 = ps(1,2);
xmax = x1;
ymax = y1;
xned = NaN;
for i=ps(2:end,:)'
    x2 = i(1);
    y2= i(2);
    if (ymax < y2) % has started to decrease
        xmax = x2;
        ymax = y2;
    end
    if (y2 < 0 && y1 > 0)
        xned = x2-y2*(x2-x1)/(y2-y1);
    end
    x1=x2;
    y1=y2;
end

fprintf("[xmax, ymax] = [%f, %f]\nxned = %f\n", xmax, ymax, xned);


%% Kvadratisk interpolation

x1 = ps(1,1);
y1 = ps(1,2);
xmax = x1;
ymax = y1;
xned = NaN;
for i=2:2:length(x)
    x2 = ps(i,1);
    y2 = ps(i,2);
    x3 = ps(i+1,1);
    y3 = ps(i+1,2);

    xv = [x1,x2,x3];
    yv = [y1,y2,y3];

    A = [[1;1;1],xv', xv'.^2];

    z = A'*A \ A'*yv';
    a=z(1);
    b=z(2);
    c=z(3);

    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x,y,vx,vy]=kastbana(h)

%KASTBANA(H) beräknar banan för ett kast med en liten boll.
%
%   Dynamiken ges av en ODE som inkluderar effekten av luftmotståndet,
%
%      r'' = -g*ez-sigma*r'*|r'|/m.
%
%   Funktionen beräknar bollens position och hastighet vid
%   tidpunkter separerade med en given steglängd. Bollen kastas från
%   (X,Y)=(0,2) med hastigheten 30 m/s i 45 graders vinkel uppåt.
%
%   Syntax:
%
%   [T,X,Y,VX,VY] = KASTBANA(H)
%
%   H       - Steglängden mellan tidpunkterna.
%   T       - Vektor med tidpunkter där bollens position och hastighet beräknats.
%   X, Y    - Vektorer med bollens x- och y-koordinater vid tidpunkterna.
%   VX, VY  - Vektorer med bollens hastigheter i x- och y-led vid tidpunkterna.

%% Tennisboll, specifikationer

m = 56e-3;     % Massan (kg) = 56 gram
ra = 6.6e-2/2; % 6.6 cm in diameter

g=9.81;      % Tyngdaccelerationen (m/s^2)

rho=1.2;     % Luftens densitet (kg/m^3)
A=ra^2*pi;   % Kroppens tvärsnittsarea (m^2)
Cd=0.47;     % Luftmotståndskoefficient,
% "drag coefficient" (dimensionslös)
% Läs mer på http://en.wikipedia.org/wiki/Drag_coefficient

sigma = rho*A*Cd/2; % Totala luftmotståndet

T  = 5;      % Sluttid
v0 = 32;     % Utkasthastighet
al = pi/4;   % Utkastvinkel

% Begynnelsevärden

r0 = [0 2]';                   % Position
r1 = [v0*cos(al) v0*sin(al)]'; % Hastighet

% ODEns högerled

f = @(u) [u(3:4); -u(3:4)*norm(u(3:4),2)*sigma/m - [0;g]];  % RHS

u = [r0;r1];
U = u';
t = 0:h:T;

% Runge-Kutta 4

for tn=t(1:end-1)
    s1 = f(u);
    s2 = f(u + h/2*s1);
    s3 = f(u + h/2*s2);
    s4 = f(u + h*s3);
    u = u + h/6*(s1 + 2*s2 + 2*s3 + s4);
    U = [U; u'];
end

x  = U(:,1);
y  = U(:,2);
vx = U(:,3);
vy = U(:,4);

end

