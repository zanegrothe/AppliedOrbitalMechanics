% Zane Grothe
% AERO 6330
% HW 3 #1
% 2/28/22

% Converting orbital elements into position(r) and velocity(v) in 3
% dimensions.
% First, convert scalar elements to perifocal coordinate system
% Second, rotate to inertial reference frame


function [rf,vf]=COE2RV(a,mu,p,ec,nu,w,in,O)

r=p/(1+ec*cos(nu)); % r scalar

% r and v vectors in perifocal coordinate system
rPQW=r*[          cos(nu) ;                 sin(nu) ; 0];
vPQW=[-sqrt(mu/p)*sin(nu) ; sqrt(mu/p)*(ec+cos(nu)) ; 0];

% Rotational matrices for inertial reference frame
Rw  = [cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1]; % ROT 3 (w)
if in < pi
    Rin = [1,0,0;0,cos(in),sin(in);0,-sin(in),cos(in)]; % ROT 1 (in)
else
    Rin = [1,0,0;0,cos(in),-sin(in);0,sin(in),cos(in)]; % ROT 1 (in)
end
RO  = [cos(O),sin(O),0;-sin(O),cos(O),0;0,0,1]; % ROT 3 (O)

A=Rw*Rin*RO; % Combined

% Rotated
rIJK=A.'*rPQW;
rf=rIJK;

vIJK=A.'*vPQW;
vf=vIJK;

