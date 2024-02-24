% Zane Grothe
% AERO 6330
% HW 4
% 3/15/22

% Function file for finding the position of L1 Lagrange Point
% Input:
%           mu   = mass parameter
%           tol  = Newton's Method tolerance limit
% Output:
%           x0   = position of L1

function x0 = L1_Position(mu,tol)
syms x
% For all solutions (5 of them):
% dUdx = x - (1-mu)*(x+mu)/sqrt((x+mu)^2)^3 - mu*(x+mu-1)/sqrt((x+mu-1)^2)^3
% For position of L1 x must be >-mu and <1-mu so:
dUdx = x - (1-mu)/(x+mu)^2 + mu/(x+mu-1)^2;
d2Udx2 = diff(dUdx); % Derivative of position of L1

fraction = 1;
x0=0.5; % Intial value for position

while abs(fraction) > tol
    
    fraction = double(subs(dUdx,x,x0))/double(subs(d2Udx2,x,x0)); % Update equation
    x0 = x0 - fraction; % Newton's Method
    
end
end

