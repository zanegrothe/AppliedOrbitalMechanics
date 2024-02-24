% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

clear all
close all
clc

% Problem 1 ~~~~~~~~~~~~~~~~~~~~
% Part b)
format long


M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass
JC=3.18; % Selected Jacobi Constant in specified range

% Nondimensional mass parameters
mu=M2/M;

% Initial Conditions
disp('Initial Conditions:')
x=-0.5;
y=0;
vx=0.566;
vy=-.904;
disp(sprintf('x = %.3f',x))
disp(sprintf('y = %.3f',y))
disp(sprintf('vx = %.3f',vx))
disp(sprintf('vy = %.3f',vy))

x0=[x,y,vx,vy];

% Numerically integrate
%tspan=[0 10]; % To see it stay with Moon
tspan=[0 50]; % To see it go back to Earth
options=odeset('RelTol',1e-6,'AbsTol',1e-8); % Set tolerences
[t,x]=ode45('CR3BP_EOM',tspan,x0,options,mu);

% Plots
colormap(cool);

% Zero velocity curves
figure(1)
[X,Y]=meshgrid(-1.5:0.01:1.5);

d1=((X + mu).^2+Y.^2).^(1/2);
r1=((X + mu-1).^2+Y.^2).^(1/2);
U=2*((1/2)*(X.^2 + Y.^2)+(1-mu)./d1+mu./r1);

contour(X,Y,U,[JC,JC],'r')
hold on
title('Earth-Moon Trajectory (nondimensional)')
xlabel('X')
ylabel('Y')
axis equal
axis square

% Planets
plot(-mu,0,'ko','MarkerSize',8,'MarkerFaceColor','g')
plot(1-mu,0,'ko','MarkerSize',3,'MarkerFaceColor','b')

% Trajectory
plot(x(:,1), x(:,2), 'k')


figure(2)
d2=sqrt((mu+x(:,1)).^2+(x(:,2)).^2);
r2=sqrt((x(:,1)-(1-mu)).^2+(x(:,2)).^2);

% Compute Jacobi Constant
C=-(x(:,3).^2+x(:,4).^2)+((x(:,1).^2+x(:,2).^2)+2*(1-mu)./d2+2*mu./r2);

plot(t,C)
ylim([3,4])
title(sprintf('Jacobi Constant = %.2f',JC))
xlabel('Nondimensional Time')
ylabel('Jacobi Constant')

EM=384467; % Distance between Earth and Moon (km)
RE=6378.1; % Radius of Earth
RM=1737.4; % Radius of Moon

% Minimum distance from Earth
EP=min(d2)*EM-RE;
disp(sprintf('Earth Periapsis = %.f km',EP))

% Minimum distance from Moon
MP=min(r2)*EM-RM;
disp(sprintf('Moon Periapsis = %.f km',MP))

