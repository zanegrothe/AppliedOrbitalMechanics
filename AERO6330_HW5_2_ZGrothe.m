% Zane Grothe
% AERO 6330
% HW 5
% 3/29/22

clear all
close all
clc
format long

% Problem 2 ~~~~~~~~~~~~~~~~~~~~


%% GIVENS

M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass
tol=1e-8; % Tolerance

% Nondimensional mass parameters
mu=M2/M;


%% FIND L1
L1x=L1_Position(mu,tol); % L1x = scalar
L1y=0;


%% FIND Trajectory IC's
syms x y

d=sqrt((x+mu)^2+y^2); % Distance from Earth to Spacecraft
r=sqrt((x+mu-1)^2+y^2); % Distance from Moon to Spacecraft

U=(x^2+y^2)/2+(1-mu)/d+mu/r; % Pseudo-potential

Uxx=diff(U,x,2); % Second derivatives
Uyy=diff(U,y,2);

Uxx=double(subs(Uxx,[x,y],[L1x,0])); % Solve second derivatives
Uyy=double(subs(Uyy,[x,y],[L1x,0]));

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
Linear_Solution_Initial_Velocity = et0_dot
xy0=[L1x+ep0,et0,ep0_dot,et0_dot];

% Period
P=2*pi/s;

% Numerically integrate
tspan=[0 P];
options=odeset('RelTol',1e-6,'AbsTol',1e-8); % Set tolerences
[t,xy]=ode45('CR3BP_EOM',tspan,xy0,options,mu);


%% PLOTS

figure(1) % (follow trajectory no matter how big)
% Planets
plot(-mu,0,'ko','MarkerSize',8,'MarkerFaceColor','g') % Green Earth
hold on
plot(1-mu,0,'ko','MarkerSize',3,'MarkerFaceColor','m') % Magenta Moon
% L1
plot(L1x,0,'s','MarkerSize',8,'color','k') % L1 point on black box
% Trajectory
plot(xy(:,1), xy(:,2), 'r')
title(sprintf('Full View Nonlinear Dynamics for %.2f from L1',delta))
xlabel('X')
ylabel('Y')
axis square
axis equal
hold off

figure(2) % (now zoom in on L1)
% L1
plot(L1x,0,'s','MarkerSize',8,'color','k') % L1 point on black box
hold on
% Trajectory
plot(xy(:,1), xy(:,2), 'r')
xlim([L1x*(1-5*abs(delta)/L1x),L1x*(1+5*abs(delta)/L1x)])
ylim([et0-L1x*5*(abs(delta)/L1x),et0+L1x*5*(abs(delta)/L1x)])
title(sprintf('Scaled View Nonlinear Dynamics for %.2f from L1',delta))
xlabel('X')
ylabel('Y')
axis square
hold off

figure(3) % Are x and y periodic?
plot(xy(:,1),'r')
hold on
plot(xy(:,2),'b')
title(sprintf('Range of Position Movements for %.2f from L1',delta))
xlabel('Nondim. Time')
ylabel('Position Change')
legend({'x','y'},'Location','east')

