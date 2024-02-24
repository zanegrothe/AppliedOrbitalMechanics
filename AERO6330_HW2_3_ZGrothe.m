% Zane Grothe
% AERO 6330
% HW 2
% 2/14/22

clear all
close all
clc

% Problem 3 ~~~~~~~~~~~~~~~~~~~~

j2=-1.08263*10^-3; % J2 dimensionless coefficient
j3=2.53244*10^-6; % J3 dimensionless coefficient
RE=6378; % Radius of Earth (km)
mu=398600; % Earth's Gravitational Parameter (km^3/s^2)

%J2=-j2*RE^2*mu; % J2 (km^5/s^2)
%J3=-j3*RE^3*mu; % J3 (km^6/s^2)

alt=100; % Altitude (km)
a=RE+alt; % Semi-major axis (km)
n=sqrt(mu/a^3); % Mean motion (/s)
in=90*pi/180; % Inclination (radians)
w=90*pi/180; % Argument of periapsis (radians)

% Calculate ----------

syms ec J2 J3 p t

A=(-3*n*J3*RE^3*sin(in))/(2*a^3*(1-ec^2)^2);
B=1-5*sin(in)^2/4;

eq1=A*B*cos(w);

%C=(3*n*J2*RE^2)/(a^2*(1-ec^2)^2);
C=(3*n*J2*RE^2)/(p^2);
%D=(J3*RE)/(2*J2*a*(1-ec^2));
D=(J3*RE)/(2*J2*p);
E=(sin(t)^2-ec*cos(t)^2)/sin(t);
F=sin(w)/ec;

eq2=C*B*(1+D*E*F);

%eqns=[eq1==0,eq2==0];
%vars=[w ec];

%[solom,solec]=solve(eqns,vars)
solve(eq2,ec)

