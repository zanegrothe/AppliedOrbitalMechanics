% Zane Grothe
% AERO 6330
% HW 5
% 3/29/22

% Function file for correcting linear solution initial velocity
% INPUT
%       r0       = initial position x;y
%       v0       = initial velocity vx;vy
%       tf       = final time before corrections
%       mu       = mass parameter
%       max_it   = maximum iterations allowed
%       toll     = tollerance of corrector
% OUTPUT
%       Xd0star  = corrected initial design variables 2x1
%       err      = position err
%       n_it     = number of iterations

function [Xd0star,err,n_it]=targetor_corrector(r0,v0,tf,mu,max_it,tol)

% Set initial values
n_it = 0;
err = 1e14;
corr = [0;0];

% Design variables
Xd0 = [v0(2);tf];

% Need error less than tollerance and iteration count below max allowed
while norm(err) > tol && n_it <= max_it
    n_it = n_it + 1; % Count it
    Xd0 = Xd0 + corr; % Make corrections
    
    x0 = [r0',v0(1),Xd0(1),reshape(eye(4,4),1,16)]; % IC array
    options=odeset('RelTol',1e-12,'AbsTol',1e-14); % Set tolerences
    
    % Calculate State Transition Matrix
    tspan=[0 tf];
    [TOUT,XOUT]=ode45('CR3BP_STM',tspan,x0,options,mu);
    PHI_tt0 = reshape(XOUT(5:20),4,4);
    PHI_24 = PHI_tt0(2,4);
    PHI_34 = PHI_tt0(3,4);
    
    % Numerically integrate new ICs
    dXOUT = CR3BP_EOM(Xd0(2),XOUT(end,1:4),options,mu);
    
    % Vy0 correction
    corr(1) = -XOUT(end,3)/(PHI_34-PHI_24*dXOUT(3)/XOUT(end,4));
    % Time correction
    corr(2) = -PHI_24/XOUT(end,4)*corr(1);
    % Error
    err = [-XOUT(end,2);-XOUT(end,3)];
end
% Output
Xd0star = Xd0;

end

