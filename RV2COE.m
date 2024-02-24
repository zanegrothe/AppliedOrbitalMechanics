% Zane Grothe
% AERO 6330
% HW 3 #1
% 2/28/22

% Converting position(r) and velocity(v) in 3 dimensions into orbital
% elements.
% r and v are [3x1] matrices, mu is the scalar gravitaional parameter

% If the variable name ends in a "v" it is a vector (nv="n vector")
% If it ends in an "h" it is a unit vector (nh="n hat")
% Otherwise the variable is the magnitude of that vector (n)


function [a,ec,in,O,w,nu]=RV2COE(r,v,mu)

% Define r and v vectors and magnitudes
rv=r;
r=norm(r);
vv=v;
v=norm(v);

% Define unit vectors for inertial reference frame (IJK)
ih=[1;0;0]; % "i hat"
jh=[0;1;0]; % "j hat"
kh=[0;0;1]; % "k hat"

% Angular momentum
hv=cross(rv,vv);
h=norm(hv);

% Node
nv=cross(kh,hv);
n=norm(nv);
nh=nv/n;

% Eccentricity (from definition)
ecv=((v^2-mu/r)*rv-dot(rv,vv)*vv)/mu;
ec=norm(ecv);
% Floating point precision inconsistancy correction
if ec < 10^-10
    ecv=zeros(3,1);
    ec=norm(ecv);
end
if abs(1-ec) < 10^-10
    ec=1;
end

% Semi-Major Axis (from definition of energy)
z=v^2/2-mu/r; % Total energy of spacecraft
if ec == 1
    a=Inf;
else
    a=-mu/(2*z);
end

% Inclination (from definition)
in=acos(dot(hv,kh)/h);
% Quadrant check
if in > pi
    in=2*pi-in;
end

% RAAN (from definition)
O=acos(dot(nh,ih));
if dot(nh,jh) < 0
    O=2*pi-O;
end

% Argument of Periapsis (from definition)
if ec == 0
    w=0;
else
    w=acos(dot(nv,ecv)/(ec*n));
end
% Quadrant check
if dot(ecv,kh) < 0
    w=2*pi-w;
end

% True anomaly (from definition)
if ec == 0
    nu=0;
else
    nu=acos(dot(ecv,rv)/(ec*r));
end
% Quadrant check
if abs(dot(rv,vv)) < 0
    nu=2*pi-nu;
end

