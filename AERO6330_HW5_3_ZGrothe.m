% Zane Grothe
% AERO 6330
% HW 5
% 3/29/22

clear all
close all
clc
format long

% Problem 3 ~~~~~~~~~~~~~~~~~~~~


%% GIVENS

M1=5.9722*10^24; % Mass of Earth (kg)
M2=7.342*10^22; % Mass of Moon (kg)
L=384400; % Distance between Earth and Moon (km)
T=4.34210393299665; % Earth Moon period (days)
M=M1+M2; % Total mass (kg)
tol=1e-12; % Tolerance

% Nondimensional mass parameters
mu=M2/M;


%%% FIND L1
L1x=L1_Position(mu,tol); % L1x = scalar
L1y=0;


%%% FIND Trajectory IC's
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
xy0=[L1x+ep0,et0,ep0_dot,et0_dot];

Pi=2*pi/s; % Period

t=linspace(0,Pi,20); % Time span

% Linear solution
EP=ep0*cos(s*t)+(et0/B3)*sin(s*t);
ET=et0*cos(s*t)-B3*ep0*sin(s*t);

X=EP+L1x;
Y=ET+L1y;

% Numerically integrate linear solution ICs
tspan=[0 Pi];
options=odeset('RelTol',1e-12,'AbsTol',1e-14); % Set tolerences
[tu,xyu]=ode45('CR3BP_EOM',tspan,xy0,options,mu);
% Find where trajectory crosses x-axis
cross=max(find(xyu(:,2)<0))+1;
xyu=xyu(1:cross,:);
tu=tu(1:cross);
tf=tu(end);

%% Targetor / Corrector
max_it=100;
r0=xy0(1:2)';
v0=xy0(3:4)';
[Xd0star,err,n_it]=targetor_corrector(r0,v0,tf,mu,max_it,tol);

% Numerically integrate corrected ICs
%xy0c=[L1x+ep0,et0,ep0_dot,-Xd0star(1)];
xy0c=[L1x+ep0,et0,ep0_dot,-.078237];
tspan=[0 2*Xd0star(2)];
%tspan=[0 2.71];
options=odeset('RelTol',1e-12,'AbsTol',1e-14); % Set tolerences
[tc,xyc]=ode45('CR3BP_EOM',tspan,xy0c,options,mu);


%% PLOT and calculations

% L1
plot(L1x,0,'s','MarkerSize',8,'color','k') % L1 point on black box
hold on
% Trajectory
plot(X,Y,'o','color','b')
plot(xyu(:,1), xyu(:,2), 'r')
plot(xyc(:,1),xyc(:,2), 'g')
sl=7;
xlim([L1x*(1-sl*abs(delta)/L1x),L1x*(1+sl*abs(delta)/L1x)])
ylim([et0-L1x*sl*(abs(delta)/L1x),et0+L1x*sl*(abs(delta)/L1x)])
title(sprintf('Linear, Nonlinear, and Targeted Dynamics for %.2f from L1',delta))
xlabel('X')
ylabel('Y')
legend('L1','Linear','Nonlinear','Targeted')


%% PERIOD of corrected orbit
x0f=xy0c(1);
y0f=xy0c(2);
df=sqrt((x0f+mu)^2+y0f^2); % Distance from Earth to Spacecraft
rf=sqrt((x0f+mu-1)^2+y0f^2); % Distance from Moon to Spacecraft
Uxxf=1-(1-mu)/df^3+(3*(1-mu)*(x0f+mu)^2)/df^5 ...
    -mu/rf^3+(3*mu*(x0f-1+mu)^2)/(rf^5);
Uyyf=1-(1-mu)/df^3+(3*(1-mu)*y0f^2)/df^5 ...
    -mu/rf^3+(3*mu*y0f^2)/(rf^5);
B1f=2-(Uxxf+Uyyf)/2;
B2f=-Uxxf*Uyyf;
sf=sqrt(-B1f+sqrt(B1f^2+B2f));
Pf=2*pi/sf; % Period (ndim)
TOP=Pf*T; % Period (days)
disp(sprintf('The Target Orbit Period is %.3f days.',TOP))


%% ERROR
er=abs(xyc(end,1)-xyc(1,1)); % Error (ndim)
TOE=er*L; % Error (km)
disp(sprintf('The Target Orbit Error is %.3f km.',TOE))

