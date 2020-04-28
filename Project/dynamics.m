function states = dynamics(mu,bodies,dt,tf)
%DYNAMICS Simulates continuous thrust trajectory
%   mu     = central body gravitational parameter.
%   bodies = all non-central bodies in the system. 
%            the first body is the spacecraft being controlled.
%   dt     = time step.
%   tf     = total simulation time (assumes t0 = 0).

% two-body Runge-Kutta integration
states = rk4(mu,bodies,dt,tf);

end

%% Function Definitions

function [rv_dot,u] = f(rv,mu,t,uprev)
% Two-Body Gravity

N = size(rv,2);
rv_dot = zeros(size(rv));

for i = 1:N
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    rv_dot(4:6,i) = - mu * rv(1:3,i) / norm(rv(1:3,i))^3;
    
    % continuous thrust for first body
    % ----- control function goes here -----
    if i == 1
        u = controls(rv(:,i),mu,t,uprev);
        rv_dot(4:6,i) = rv_dot(4:6,i) + u;
    end
    % ----- end of control function -----
end

end

function xvect = rk4(mu,x0,dt,T)
% 4-Stage Runge-Kutta Propagator

xvect = repmat(x0,1,1,T/dt);
u = ones(3,1);

for jj = 1 : T/dt
    
    xvect(:,:,jj) = x0;
    tvect = linspace(dt*(jj-1),dt*jj,4);
    
    % ============= RK4 =============
    k1 = dt * f( x0,mu,tvect(1),u);
    k2 = dt * f((x0 + k1/2),mu,tvect(2),u);
    k3 = dt * f((x0 + k2/2),mu,tvect(3),u);
    k4 = dt * f((x0 + k3)  ,mu,tvect(4),u);
    x0 = x0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    [~,u] = f( x0,mu,tvect(1),u);
    
end

end