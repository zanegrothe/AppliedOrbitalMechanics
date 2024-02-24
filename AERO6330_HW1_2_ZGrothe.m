% Zane Grothe
% AERO 6330
% HW 1
% 1/27/22

clear all
close all
clc

rE = 6378; % Radius of Earth (km)
a = 7*rE; % Semi-Major Axis (km)
e = 0.92; % Eccentricity
tol = 10^-12; % Tolerance (radians)

M = 2*pi/10; % Initial value for Mean Anomaly (radians)
cvMA = 0; % Counting variable for Mean Anomaly loop

while M <= 2*pi
    
    syms E
    K = E - e*sin(E) - M; % Kepler's Equation
    dK = diff(K); % Derivative of Kepler's Equation
    
    fraction = 1;
    cvNM = 0; % Counting variable for Newton's Method loop
    E0 = M; % Initial value for Eccentric Anomaly
    
    while abs(fraction) > tol
        
        fraction = double(subs(K,E,E0))/double(subs(dK,E,E0));
        E0 = E0 - fraction; % Newton's Method
        cvNM = cvNM + 1; % Increase counting variable
        
    end
    
    M = M*180/pi; % Convert Mean Anomaly from radians to degrees
    E0 = E0*180/pi; % Convert Eccentric Anomaly from radians to degrees
    I = cvNM; % Total number of iterations for convergence
    Difference = abs(M-E0); % Difference between Anomalies (degrees)
    
    if cvMA == 0;
        % Initial matrix values for...
        M_m = M; % Mean Anomaly
        E_m = E0; % Eccentric Anomaly
        I_m = I; % Iterations
        D_m = Difference; % Difference
    else
        % Additional matrix values for...
        M_m = [M_m,M]; % Mean Anomalies
        E_m = [E_m,E0]; % Eccentric Anomalies
        I_m = [I_m,I]; % Iterations
        D_m = [D_m,Difference]; % Differences
    end
    
    M = M*pi/180; % Convert Mean Anomaly back to radians
    M = M + 2*pi/10; % Increase Mean Anomaly by one tenth of orbit
    cvMA = cvMA + 1; % Increase counting variable
    
end

Data = [M_m',E_m',I_m',D_m']; % Arrange outputs into matrix

disp('  Mean An.  Ecc. An.    Iter.    A. Diff.')
disp(Data)

