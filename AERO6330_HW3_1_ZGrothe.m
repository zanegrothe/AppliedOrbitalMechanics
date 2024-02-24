% Zane Grothe
% AERO 6330
% HW 3
% 2/28/22

clear all
close all
clc

% Problem 1 ~~~~~~~~~~~~~~~~~~~~

% Implement a function COE2RV to convert from Keplerian orbit elements to
% position and velocity vector in planet centered inertial coordinates.
% Implement a function RV2COE to perform the inverse operation (from
% position/velocity vector in planet centered inertial coordinates to
% Keplerian orbit elements).


ec=; % Eccentricity
% Semi-Major Axis or Periapsis(km)
SP=; % (km)
if ec == 1
    rp=SP; % Periapsis for a parabolic orbit (km)
else
    a=SP; % Semi-major axis (km)
end
in=; % Inclination (deg)
O=; % RAAN (deg)
w=; % Argument of periapsis (deg)
nu=; % True anomaly (deg)

mu=398600; % Gravitational mass parameter (km^3/s^2)

% Convert angles radians
in=in*pi/180;
O=O*pi/180;
w=w*pi/180;
nu=nu*pi/180;

% Check if true anomaly is less than asymptote angle
if ec > 1
    if nu > acos(-1/ec) || nu < -acos(-1/ec)
        disp('Invalid true anomaly!')
        disp('You have a hyperbolic trajectory (eccentricity is greater than 1)')
        disp('A valid true anomaly for this trajectory must be less than:')
        disp(acos(-1/ec)*180/pi)
        disp('and greater than')
        disp(-acos(-1/ec)*180/pi)
        return
    end
end

% Parameter (km) and display movement identification
if ec == 0
    p=a*(1-ec^2);
    disp('Position(km) and Velocity(km/s) vectors for your circular orbit:')
elseif 0 < ec && ec < 1
    p=a*(1-ec^2);
    disp('Position(km) and Velocity(km/s) vectors for your elliptical orbit:')
elseif ec == 1
    p=2*rp;
    disp('Position(km) and Velocity(km/s) vectors for your parabolic trajectory:')
else
    p=a*(ec^2-1);
    disp('Position(km) and Velocity(km/s) vectors for your hyperbolic trajectory:')
end 

% Convert Orbital Elements to r and v
if ec == 1
    [r,v]=COE2RV(rp,mu,p,ec,nu,w,in,O)
else
    [r,v]=COE2RV(a,mu,p,ec,nu,w,in,O)
end

% Covert r and v back to Orbital Elements
% If the variable ends in "t" it is temporary
[a,ec,in,O,wt,nut]=RV2COE(r,v,mu);

% Display movement identification
if ec == 0
    disp('Orbital Elements for your circular orbit:')
elseif 0 < ec && ec < 1
    disp('Orbital Elements for your elliptical orbit:')
elseif ec == 1
    disp('Orbital Elements for your parabolic trajectory:')
else
    disp('Orbital Elements for your hyperbolic trajectory:')
end 

a
ec

% Defining arbitraty (originally specified) elements for a circular orbit
if ec == 0
    w=wt+w;
    nu=nut+nu;
end

% Convert angles back to degrees
in=in*180/pi
O=O*180/pi
w=w*180/pi
nu=nu*180/pi

