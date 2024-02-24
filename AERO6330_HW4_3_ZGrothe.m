% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

clear all
close all
clc

% Problem 3 ~~~~~~~~~~~~~~~~~~~~
format short

mu=0.0121; % Mass parameter

% Initial conditions (L1 found from previous problem)
x0=0.8369;
y0=0;

% Find roots of characteristic equation
roots=L1_Lyapunov(mu,x0,y0)

% Check for stability
stability=isreal(roots(1,1));
if stability == 1
    disp('A real root has been found. This position is unstable.')
else
    disp('A real root has not been found. This position is stable.')
end

