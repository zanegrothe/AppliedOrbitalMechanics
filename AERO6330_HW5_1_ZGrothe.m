% Zane Grothe
% AERO 6330
% HW 5
% 3/29/22

clear all
close all
clc
format long

% Problem 1 ~~~~~~~~~~~~~~~~~~~~


%% GIVENS

G=1; % Nondimensional Gravitational Constant
M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass
tol=1e-8; % Tolerance

% Nondimensional mass parameters
mu=G*M2/M;


%% FIND L1
L1x=L1_Position(mu,tol); % L1x = scalar
L1y=0;


%% FIND TRAJECTORY IC's
syms x y

d=sqrt((x+mu)^2+y^2); % Distance from Earth to Spacecraft
r=sqrt((x+mu-1)^2+y^2); % Distance from Moon to Spacecraft

U=(x^2+y^2)/2+(1-mu)/d+mu/r; % Pseudo-potential

Uxx=diff(U,x,2); % Second derivatives
Uyy=diff(U,y,2);

Uxx=double(subs(Uxx,[x,y],[L1x,0])); % Solve second derivatives
Uyy=double(subs(Uyy,[x,y],[L1x,0]));
del=sqrt(Uyy/Uxx);

% Constants
B1=2-(Uxx+Uyy)/2;
B2=-Uxx*Uyy;
s=sqrt(B1+sqrt(B1^2+B2));
B3=(s^2+Uxx)/2/s;

% L1 offset
delta=0.01;

% Initial Conditions
ep0      = delta;
et0      = 0;
ep0_dot  = et0*s/B3;
et0_dot  = -ep0*s*B3;
xy0=[ep0,et0,ep0_dot,et0_dot];


%% FIND Roots of Lyapunov Stability Analysis
roots=L1_Lyapunov(mu,L1x,L1y); % roots = [lam1;lam2;lam3;lam4];


%% FIND Linear Solution Trajectory
alpha1=(roots(1,1)^2-Uxx)/2/roots(1,1);
alpha2=(roots(2,1)^2-Uxx)/2/roots(2,1);
alpha3=(roots(3,1)^2-Uxx)/2/roots(3,1);
alpha4=(roots(4,1)^2-Uxx)/2/roots(4,1);

t0=0; % Initial time

A1=exp(-roots(1,1)*t0)/(roots(1,1)^2-roots(3,1)^2)...
    * (-ep0*alpha3*roots(3,1) - ep0_dot*alpha3*del...
    + et0*roots(3,1)*del + et0_dot);
A2=exp(roots(1,1)*t0)/(roots(1,1)^2-roots(3,1)^2)...
    * (-ep0*alpha3*roots(3,1) + ep0_dot*alpha3*del...
    - et0*roots(3,1)*del + et0_dot);
A3=exp(-roots(3,1)*t0)/(roots(1,1)^2-roots(3,1)^2)...
    * (ep0*alpha1*roots(1,1) + ep0_dot*alpha1*del...
    - et0*roots(1,1)*del - et0_dot);
A4=exp(roots(3,1)*t0)/(roots(1,1)^2-roots(3,1)^2)...
    * (ep0*alpha1*roots(1,1) - ep0_dot*alpha1*del...
    + et0*roots(1,1)*del - et0_dot);

P=2*pi/s; % Period

t=linspace(0,P); % Time span

% Original
%EP=A3*exp(roots(3,1)*t) + A4*exp(roots(4,1)*t);
%ET=A3*alpha3*exp(roots(3,1)*t) - A4*alpha4*exp(roots(4,1)*t);
% Transformed
EP=ep0*cos(s*t)+(et0/B3)*sin(s*t);
ET=et0*cos(s*t)-B3*ep0*sin(s*t);

% Coordinate frame adjustment
X=EP+L1x;
Y=ET+L1y;

%% PLOTS

figure(1) % Follow trajectory no matter how big
% Planets
plot(-mu,0,'ko','MarkerSize',8,'MarkerFaceColor','g') % Green Earth
hold on
plot(1-mu,0,'ko','MarkerSize',3,'MarkerFaceColor','m') % Magenta Moon
% L1
plot(L1x,0,'s','MarkerSize',8,'color','k') % L1 point on black box
% Trajectory
plot(X,Y,'b')
title(sprintf('Full View Linear Dynamics for %.2f from L1',delta))
xlabel('X')
ylabel('Y')
axis square
axis equal
hold off

figure(2) % Zoomed in on L1
% L1
plot(L1x,0,'s','MarkerSize',8,'color','k') % L1 point on black box
hold on
% Trajectory
plot(X,Y,'b')
xlim([L1x*(1-5*abs(delta)/L1x),L1x*(1+5*abs(delta)/L1x)])
ylim([et0-L1x*5*(abs(delta)/L1x),et0+L1x*5*(abs(delta)/L1x)])
title(sprintf('Scaled View Linear Dynamics for %.2f from L1',delta))
xlabel('X')
ylabel('Y')
axis square
hold off

figure(3) % Are Epsilon and Eta periodic?
plot(EP,'r')
hold on
plot(ET,'b')
title(sprintf('Range of Position Movements for %.2f from L1',delta))
xlabel('Nondim. Time')
ylabel('Position Change')
legend({'Epsilon','Eta'},'Location','east')

