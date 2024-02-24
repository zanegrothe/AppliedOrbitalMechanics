% Zane Grothe
% AERO 6330
% HW 5
% 3/29/22

% Function file for finding the state transition matrix for cr3bp
% INPUT
%       position vector (x,y):        xx(1),xx(2)
%       velocity vector (vx,vy):      xx(3),xx(4)
%       mass parameter:               mu
% OUTPUT
%       velocity vector (vx,vy):      yy(1),yy(2)
%       acceleration vector (ax,ay):  yy(3),yy(4)

function yy=CR3BP_STM(t,xx,options,mu)
% Position
x = xx(1);
y = xx(2);
% Velocity
vx = xx(3);
vy = xx(4);
% Derivatives
dx = vx;
dy = vy;
% Distances
d = sqrt((x+mu).^2 + y.^2);
r = sqrt((x+mu-1).^2 + y.^2);
% Accelerations
ax=x+2*vy+(1-mu).*(-mu-x)./(d.^3)+mu.*(1-mu-x)./(r.^3); 
ay=y-2*vx-(1-mu).*y./(d.^3)-mu.*y./(r.^3);

% Creating output vector
yy = zeros(20,1);
yy(1) = dx;
yy(2) = dy;
yy(3) = ax;
yy(4) = ay;

% State Transition Matrix format
%[1 5 9  13
% 2 6 10 14
% 3 7 11 15
% 4 8 12 16]

% Create PHI matrix
PHI_tt0 = reshape(xx(5:20),4,4);

OMEGA = [0,2;-2,0];

Uxx=1-(1-mu)./d.^3+(3*(1-mu)*(x+mu).^2)./d.^5 ...
    -mu/r.^3+(3*mu*(x-1+mu).^2)./(r.^5);
Uxy=3*(1-mu)*(x+mu).*y./d.^5+3*mu*(x-1+mu).*y./r.^5;
Uyy=1-(1-mu)./d.^3+(3*(1-mu)*y.^2)./d.^5 ...
    -mu/r.^3+(3*mu*y.^2)./(r.^5);
Uyx=Uxy;

U = [Uxx,Uxy;Uyx,Uyy];

% Adjustment and matrix
A = [zeros(2),eye(2);U,OMEGA];
dPHI_tt0 = A*PHI_tt0;

% Reformat output vector
yy(5:20) = reshape(dPHI_tt0,16,1);
end

