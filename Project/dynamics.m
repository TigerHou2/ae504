function [states,dv] = dynamics(mu,bodies,dt,t0,tf,uvect,tvect)
%DYNAMICS Simulates continuous thrust trajectory
%   mu     = central body gravitational parameter.
%   bodies = all non-central bodies in the system. 
%            the first body is the spacecraft being controlled.
%   dt     = time step.
%   tf     = total simulation time (assumes t0 = 0).

% two-body Runge-Kutta integration
[states,dv] = rk4(mu,bodies,dt,t0,tf,uvect,tvect);

end

%% Function Definitions

function [rv_dot,u] = f(rv,mu,t,uvect,tvect)
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
        u = controls(rv(1:3,i),mu,uvect,tvect,t);
        rv_dot(4:6,i) = rv_dot(4:6,i) + u;
    end
    % ----- end of control function -----
end

end

function [xvect,dv] = rk4(mu,x0,dt,t0,tf,uu,tt)
% 4-Stage Runge-Kutta Propagator

xvect = repmat(x0,1,1,(tf-t0)/dt);
tvect = linspace(t0,tf,(tf-t0)/dt+1);
dv = 0;

for jj = 1 : (tf-t0)/dt
    
    xvect(:,:,jj) = x0;
    ttvect = linspace(tvect(jj),tvect(jj+1),4);
    u = controls(x0(1:3),mu,uu,tt,ttvect(1));
    dv = dv + norm(u)*dt;
    
    % ============= RK4 =============
    k1 = dt * f( x0,mu,ttvect(1),uu,tt);
    k2 = dt * f((x0 + k1/2),mu,ttvect(2),uu,tt);
    k3 = dt * f((x0 + k2/2),mu,ttvect(3),uu,tt);
    k4 = dt * f((x0 + k3)  ,mu,ttvect(4),uu,tt);
    x0 = x0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end

end