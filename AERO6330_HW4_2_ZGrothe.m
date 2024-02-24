% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

clear all
close all
clc

% Problem 2 ~~~~~~~~~~~~~~~~~~~~

tol = 10^-8; % Tolerance (radians)

M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass

mu=M2/M; % Nondimensional mass parameters


% Find postion of L1
x0=L1_Position(mu,tol);


disp(sprintf('Nondimensional position of L1 is %.4f',x0))
EM=384467; % Distance between Earth and Moon (km)
L1=x0*EM;
disp(sprintf('Position of L1 is %.f km right of barycenter',L1))


d1=x0+mu;
r1=1-x0-mu;

% Compute the Jacobi Constant
C=x0^2+2*(1-mu)/d1+2*mu/r1;
disp(sprintf('The Jacobi Constant is %.3f',C))

M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass

[xm,ym]=meshgrid(-1.5:.01:1.5);

x1=-mu; % x coordinate of Earth
x2=1-mu; % x coordinate of Moon
d=sqrt((xm-x1).^2+ym.^2); % Distance from Earth to Spacecraft
r=sqrt((xm-x2).^2+ym.^2); % Distance from Moon to Spacecraft
U=-(xm.^2+ym.^2)/2-(1-mu)./d-mu./r; % Pseudo-potential

% Plot
p=C/-2; % Peak level

colormap(cool);
contourf(xm,ym,U,[p p])
hold on
plot(-mu,0,'ko','MarkerSize',8,'MarkerFaceColor','g')
plot(1-mu,0,'ko','MarkerSize',3,'MarkerFaceColor','b')
plot(x0,0,'x','MarkerSize',8,'color','r') % L1 point on red x
title(sprintf('Jacobi Constant of %.3f',C))
axis equal
axis square



