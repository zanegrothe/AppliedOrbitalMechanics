% Zane Grothe
% AERO 6330
% HW 6
% 4/19/22

clear all
close all
clc

%% Problem 1 ~~~~~~~~~~~~~~~~~~~~

% Givens

% Colored arcs converted to radians
red=35*pi/180;      % red
green=45*pi/180;    % green
blue=30*pi/180;     % blue

% Associated angles
RED=red;    % Red
BLUE=blue;  % Blue

gp=pi/2-green; % g' arc is 90 degrees minus green arc

rp=asin(sin(RED)*sin(gp)/sin(BLUE)); % r' arc is found from sine rule

% c arc is found from Napier's analogies
c=2*atan(tan(.5*(gp+rp))*cos(.5*(BLUE+RED))/cos(.5*(BLUE-RED)));

purple=pi/2-c;          % purple arc is 90 degrees minus c arc
purple=purple*180/pi;   % convert back to degrees

% Display results
disp('Problem 1:')
disp(sprintf('The latitude of our field agent is %.3f degrees North.',purple))
disp(' ')


%% Problem 2 ~~~~~~~~~~~~~~~~~~~~

% Givens

lam=10*pi/180;          % Target central angle (10deg)
phiE=2*pi-80*pi/180;    % Target azimuth angle (80deg) reference North COUNTERCLOCKWISE
% The equation for azimuth is based on a CLOCKWISE measured angle so we
% need 2pi-"azimuth" angle to calculate

LAT_T=32*pi/180;    % Latitude of Target (32deg North)
LON_T=85*pi/180;    % Longitude of Target (85deg West)
LAT_pT=pi/2-LAT_T;  % Latitude' of Target (from North pole)

% Calculations

% Solve for Latitude' of SSP (from North pole)
syms latpssp
phi_E=acos((cos(LAT_pT)-cos(latpssp)*cos(lam))/(-sin(latpssp)*sin(lam)))+pi==phiE;
f=solve(phi_E,latpssp);
LAT_pSSP=real(double(f(2,1)));

LAT_SSP=pi/2-LAT_pSSP;  % Latitude of SSP
LAT_SSP=LAT_SSP*180/pi; % convert back to degrees

% Solve for angle between Target and SSP
deltaL=acos((cos(lam)-cos(LAT_pT)*cos(LAT_pSSP))/(sin(LAT_pT)*sin(LAT_pSSP)));

LON_SSP=LON_T-deltaL;   % Subtract latitude difference from target latitude
LON_SSP=LON_SSP*180/pi; % convert back to degrees

% Display results
disp('Problem 2 Part C):')
disp('The current location of our sub-satellite point is')
disp(sprintf('%.3f deg N,',LAT_SSP))
disp(sprintf('%.3f deg W.',LON_SSP))
disp(' ')


%% Problem 3 ~~~~~~~~~~~~~~~~~~~~

% Givens

H=200;      % Orbit altitude (km)
RE=6378.1;  % Radius of Earth (km)
mu=398600;  % Earth gravitational parameter (km^3/s^2)

% Calculations

% Find satellite's coverage
eMIN=atand(1.5/5);                  % Minumum elevation angle of target (deg)
eta=asind(RE*sind(eMIN+90)/(RE+H)); % Nadir angle of target (deg)
lamMAXa=90-eMIN-eta;                % Maximum central angle of target (deg)
lamMAXd=RE*lamMAXa*2*pi/180;        % Arc diameter of access area (km)

% Find relative speeds and distances covered by target and satellite
RExT=RE*cos(LAT_T); % Target's distance from Earth's rotational axis (km)
vT=2*pi*RExT/24;    % Target's tangental velocity (km/hr)

vSo=sqrt(mu/(RE+H))*3600;   % Satellite's tangental velocity in orbit (km/hr)
wS=vSo/(RE+H);              % Satellite's angular velocity (rad/hr)
vS=wS*RE;                   % Satellite's tangental velocity on Earth (km/hr)

SATvh=sqrt(vT^2+vS^2);  % Target speed relative to satellite (km/hr)
SATah=atand(vS/vT);     % Target direction relative to satellite (deg)

% Find orbit planes
xTd=lamMAXd*cosd(SATah);    % Distance between orbit planes (km)
xTa=xTd/RE*180/pi;          % Angle between orbit planes (deg)

yTd=lamMAXd*sind(SATah);    % Distance between satellite true anomalies (km)
yTa=yTd/RE*180/pi;          % Angle between satellite true anomalies (deg)

P=360/xTa; % Number of orbit planes

% Display results
disp('Problem 3:')
disp(sprintf('We will need %.f orbit planes for continuous coverage (just for one day).',P))
disp(sprintf('The angle between each consecutive orbit plane is %.4f degrees.',xTa))
disp(sprintf('The angle between each consecutive satellite latitude is %.4f degrees.',yTa))


% Propagate orbit (equally spaced orbits) ~~~~~~~~~~~~~~~~~~~~
ec=0;                   % Eccentricity
a=RE+H;                 % Semi-major axis (km)
in=90;                  % Inclination (deg)
O=0;                    % RAAN (deg)
dO=xTa*20;              % change in RAAN (deg)
w=32;                   % Argument of periapsis (deg)
dw=yTa*20;              % change in Argument of periapsis (deg)
num=linspace(0,360);    % True anomaly (deg)

% Convert angles radians
in=in*pi/180;
O=O*pi/180;
w=w*pi/180;
num=num*pi/180;

p=a*(1-ec^2); % Parameter

% Convert Orbital Elements to r and v
k=1;
z=1;
while O < 2*pi
    while k < 101
        nu=num(1,k);
        [r,v]=COE2RV(a,mu,p,ec,nu,w,in,O);
        if k == 1
            rm=r;
            p0=r;
        else
            rm=[rm,r]; %#ok<AGROW>
        end
        k=k+1;
    end
    if z == 1
        rtm=rm;
        p0t=p0;
    else
        rtm=[rtm;rm]; %#ok<AGROW>
        p0t=[p0t,p0]; %#ok<AGROW>
    end
    O=O+(dO*pi/180);
    w=w+(dw*pi/180);
    z=z+1;
    k=1;
end

% Plots
figure(1) % Plot equally spaced orbits
u=360/dO;
for b=1:u
    h=3*(b-1);
    plot3(rtm(1+h,:),rtm(2+h,:),rtm(3+h,:),'k')
    hold on
    if b == 1
        plot3(p0t(1,b),p0t(2,b),p0t(3,b),'ks','MarkerSize',6,'MarkerFaceColor','r')
    else
        plot3(p0t(1,b),p0t(2,b),p0t(3,b),'ks','MarkerSize',6,'MarkerFaceColor','b')
    end
end
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Z Position (km)')
title(sprintf('Orbit Trajectories for %.f Orbit Planes with Initial Satellite Positions',u))


% Propagate first 10 orbits ~~~~~~~~~~~~~~~~~~~~
ec=0;                   % Eccentricity
a=RE+H;                 % Semi-major axis (km)
in=90;                  % Inclination (deg)
O=0;                    % RAAN (deg)
dO=xTa;                 % change in RAAN (deg)
w=32;                   % Argument of periapsis (deg)
dw=yTa;                 % change in Argument of periapsis (deg)
num=linspace(0,360);    % True anomaly (deg)

% Convert angles radians
in=in*pi/180;
O=O*pi/180;
w=w*pi/180;
num=num*pi/180;

p=a*(1-ec^2); % Parameter

% Convert Orbital Elements to r and v
k=1;
z=1;
while z < P+1
    while k < 101
        nu=num(1,k);
        [r,v]=COE2RV(a,mu,p,ec,nu,w,in,O);
        if k == 1
            rm=r;
            p0=r;
        else
            rm=[rm,r]; %#ok<AGROW>
        end
        k=k+1;
    end
    if z == 1
        rtm=rm;
        p0t=p0;
    else
        rtm=[rtm;rm]; %#ok<AGROW>
        p0t=[p0t,p0]; %#ok<AGROW>
    end
    O=O+(dO*pi/180);
    w=w-(dw*pi/180);
    z=z+1;
    k=1;
end

% Plots
figure(2) % Plot first 10 orbits
for b=1:10
    h=3*(b-1);
    plot3(rtm(1+h,:),rtm(2+h,:),rtm(3+h,:),'k')
    hold on
    if b == 1
        plot3(p0t(1,b),p0t(2,b),p0t(3,b),'ks','MarkerSize',6,'MarkerFaceColor','r')
    else
        plot3(p0t(1,b),p0t(2,b),p0t(3,b),'ks','MarkerSize',6,'MarkerFaceColor','b')
    end
end
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Z Position (km)')
title(sprintf('Orbit Trajectories for First %.f Orbit Planes with Initial Satellite Positions',b))

figure(3) % Plot all satellites
for bb=1:P
    if bb == 1
        plot3(p0t(1,bb),p0t(2,bb),p0t(3,bb),'ks','MarkerSize',6,'MarkerFaceColor','r')
        hold on
    else
        plot3(p0t(1,bb),p0t(2,bb),p0t(3,bb),'ks','MarkerSize',6,'MarkerFaceColor','b')
    end
end
hold off
xlabel('X Position (km)')
ylabel('Y Position (km)')
zlabel('Z Position (km)')
title(sprintf('All Satellites for %.f Orbit Planes with Initial Satellite Positions',P))
