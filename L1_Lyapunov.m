% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

% Function file for finding the roots of the characteristic equation
% Input:
%           mu      = mass parameter
%           x0      = initial condition x
%           y0      = initial condition y
% Output:
%           roots   = roots (lambda 1-4)

function roots=L1_Lyapunov(mu,x0,y0)

syms ep et lambda

d=sqrt((ep+mu)^2+et^2); % Distance from Earth to Spacecraft
r=sqrt((ep+mu-1)^2+et^2); % Distance from Moon to Spacecraft

U=(ep^2+et^2)/2+(1-mu)/d+mu/r; % Pseudo-potential

Uxx=diff(U,ep,2); % Second derivatives
Uyy=diff(U,et,2);

% Lyapunov Stability Analysis
A=-[0  ,0  ,1 ,0;
    0  ,0  ,0 ,1;
    Uxx,0  ,0 ,2;
    0  ,Uyy,-2,0];

% Characteristic Equation (Lyapunov Stability Analysis)
s=det(lambda*eye(4,4)-A);

% Roots (together then sorted)
lm=solve(s,lambda);
lm1=lm(2,1);
lm2=lm(4,1);
lm3=lm(1,1);
lm4=lm(3,1);

% Stability of I.C.
lam1=double(subs(lm1,[ep,et],[x0,y0]));
lam2=double(subs(lm2,[ep,et],[x0,y0]));
lam3=double(subs(lm3,[ep,et],[x0,y0]));
lam4=double(subs(lm4,[ep,et],[x0,y0]));
roots=[lam1;lam2;lam3;lam4];

end

