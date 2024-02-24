% Zane Grothe
% AERO 6330
% HW 2
% 2/14/22 <3

clear all
close all
clc

% Problem 2 ~~~~~~~~~~~~~~~~~~~~

J2=1.08263*10^-3; % J2 coefficient
RE=6378; % Radius of Earth (km)
mu=398600; % Earth's Gravitational Parameter (km^3/s^2)


% Eccentricity effect on RAAN ~~~~~~~~~~

a=26554; % Semi-Major Axis (km)
ec=linspace(0.1,0.9); % Eccentricity
p=a*(1-ec.^2); % Semi-Latus Rectums (km)
n=2*pi; % Mean motion
T=2*pi*sqrt(a^3/mu); % Period (s)
S_opd=86400/T; % spacecraft orbits per day
in=asin(sqrt(4/5))+pi; % Inclination (radians)

% delta Omega i.e. procession of RAAN (radians per orbit)
DO_ec=-3/2*n*J2*RE^2*cos(in)./p.^2;

% Convert to degrees
DO_ec_deg=DO_ec*180/pi;

% Processions
E_dpd=360/365.45; % Earth degrees per day
E_dpo=E_dpd/S_opd; % Earth degrees per SC orbit

% Plot
plot(ec,DO_ec_deg,'color','b')
hold on
plot(ec,E_dpo,'color','r')
hold off
xlim([0.1,0.9])
ylim([min(DO_ec_deg),max(DO_ec_deg)])
xlabel('Eccentricity 0.1 to 0.9')
ylabel('Procession (degrees per orbit)')
title('Procession of RAAN Based on Eccentricity')
legend({'RAAN procession','Earth movement'},'Location','northeast')

