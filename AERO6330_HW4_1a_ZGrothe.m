% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

clear all
close all
clc

% Problem 1 ~~~~~~~~~~~~~~~~~~~~
% Part a)

G=1; % Nondimensional Gravitational Constant
M1=5.9722*10^24; % Mass of Earth
M2=7.342*10^22; % Mass of Moon
M=M1+M2; % Total mass
R=1; % Nondimensional distance between M1 and M2

% Nondimensional mass parameters
mu1=M1/M;
mu2=M2/M;

[x,y]=meshgrid(-2:.01:2);

omega=sqrt(mu1+mu2); % Angular velocity
x1=-mu2; % x coordinate of Earth
x2=1-mu2; % x coordinate of Moon
d=sqrt((x-x1).^2+y.^2); % Distance from Earth to Spacecraft
r=sqrt((x-x2).^2+y.^2); % Distance from Moon to Spacecraft
U=-omega^2*R/2*(x.^2+y.^2)-G*mu1./d-G*mu2./r; % Pseudo-potential

% Comb data along X axis for gateway Jacobi Constants
Pl=findpeaks(U(201,:)); % Minimums along x axis
P=[Pl(1,2),Pl(1,3),Pl(1,1)]; % Rearranging to Lagrange's order
JC=P*-2; % Compute Jacobi Constant
disp('Jacobi Constants for Earth-Moon gateway openings:')
disp('L1')
disp(JC(1,1))
disp('L2')
disp(JC(1,2))
disp('L3')
disp(JC(1,3))

[col,row]=size(P);

% Plots

% Three plots for picked Jacobi Constants
for k = 1:row
    figure(k)
    colormap(cool);
    contourf(x,y,U,[P(1,k) P(1,k)])
    hold on
    plot(-mu2,0,'ko','MarkerSize',8,'MarkerFaceColor','g') % Earth
    plot(1-mu2,0,'ko','MarkerSize',3,'MarkerFaceColor','b') % Moon
    title(sprintf('Jacobi Constant = %.4f',JC(1,k)))
    xlabel(sprintf('Gateway at L_%.0f opens',k))
    axis equal
    axis square
end

% Plot of JC in specified Range
figure(4)
colormap(cool);
contourf(x,y,U,[-1.59 -1.59])
hold on
plot(-mu2,0,'ko','MarkerSize',8,'MarkerFaceColor','g')
plot(1-mu2,0,'ko','MarkerSize',3,'MarkerFaceColor','b')
title('Jacobi Constant = 3.18')
xlabel(sprintf('Earth Moon connected at L_%.0f',1))
axis equal
axis square

