% Zane Grothe
% AERO 6330
% HW 2
% 2/14/22 <3

clear all
close all
clc

% Problem 1 ~~~~~~~~~~~~~~~~~~~~

J2=1.08263*10^-3; % J2 coefficient
RE=6378; % Radius of Earth (km)
mu=398600; % Earth's Gravitational Parameter (km^3/s^2)



% Inclination effect on RAAN ~~~~~~~~~~

a=15000; % Semi-Major Axis (km)
ec=0.8; % Eccentricity
p=a*(1-ec^2); % Semi-Latus Rectum (km)
n=2*pi; % Mean motion
T=2*pi*sqrt(a^3/mu); % Period (s)
S_opd=86400/T; % spacecraft orbits per day

in=linspace(0,2*pi); % Inclination (radians)

% delta Omega i.e. procession of RAAN (radians per orbit)
DO_in=-3/2*n*J2*RE^2*cos(in)/p^2; 

% Convert to degrees
in_deg=in*180/pi; 
DO_in_deg=DO_in*180/pi;

% Processions
E_dpd=360/365.45; % Earth degrees per day
E_dpo=E_dpd/S_opd; % Earth degrees per SC orbit

% Plots
figure(1)
plot(in_deg,DO_in_deg,'color','b')
hold on
plot(in_deg,E_dpo,'color','r')
hold off
xlim([0,360])
ylim([min(DO_in_deg),max(DO_in_deg)])
xlabel('Inclination 0 to 360 degrees')
ylabel('Procession (degrees per orbit)')
title('Procession of RAAN Based on Inclination')
legend({'RAAN procession','Earth movement'},'Location','northeast')



% Eccentricity effect on RAAN ~~~~~~~~~~

a=15000; % Semi-Major Axis (km)
ec=linspace(0.1,0.9); % Eccentricity
p=a*(1-ec.^2); % Semi-Latus Rectums (km)
n=2*pi; % Mean motion
T=2*pi*sqrt(a^3/mu); % Period (s)
S_opd=86400/T; % spacecraft orbits per day
in=150*pi/180; % Inclination (radians)

% delta Omega i.e. procession of RAAN (radians per orbit)
DO_ec=-3/2*n*J2*RE^2*cos(in)./p.^2;

% Convert to degrees
DO_ec_deg=DO_ec*180/pi;

% Processions
E_dpd=360/365.45; % Earth degrees per day
E_dpo=E_dpd/S_opd; % Earth degrees per SC orbit

% Plot
figure(2)
plot(ec,DO_ec_deg,'color','b')
hold on
plot(ec,E_dpo,'color','r')
hold off
xlim([0.1,0.9])
ylim([min(DO_ec_deg),max(DO_ec_deg)])
xlabel('Eccentricity 0 to 1')
ylabel('Procession (degrees per orbit)')
title('Procession of RAAN Based on Eccentricity')
legend({'RAAN procession','Earth movement'},'Location','northeast')



% Semi-Major Axis effect on RAAN ~~~~~~~~~~

a=linspace(15000,65000); % Perigee altitude (km)
ec=0.8; % Eccentricity
p=a*(1-ec^2); % Semi-Latus Rectum (km)
n=2*pi; % Mean motion
T=2*pi*sqrt(a.^3/mu); % Period (s)
S_opd=86400./T; % spacecraft orbits per day
in=150*pi/180; % Inclination (radians)

% delta Omega i.e. procession of RAAN (radians per orbit)
DO_a=-3/2*n*J2*RE^2*cos(in)./p.^2;

% Convert to degrees
DO_a_deg=DO_a*180/pi;

% Processions
E_dpd=360/365.45; % Earth degrees per day
E_dpo=E_dpd./S_opd; % Earth degrees per SC orbit

% Plot
figure(3)
plot(a,DO_a_deg,'color','b')
hold on
plot(a,E_dpo,'color','r')
hold off
xlim([min(a),max(a)])
ylim([min(DO_a_deg),max(DO_a_deg)])
xlabel('Semi-Major Axis (km)')
ylabel('Procession (degrees per orbit)')
title('Procession of RAAN Based on Semi-Major Axis')
legend({'RAAN procession','Earth movement'},'Location','northeast')


