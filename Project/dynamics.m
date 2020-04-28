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

function rv_dot = f(rv,mu,t)
% Two-Body Gravity

N = size(rv,2);
rv_dot = zeros(size(rv));

for i = 1:N
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    rv_dot(4:6,i) = - mu * rv(1:3,i) / norm(rv(1:3,i))^3;
    % continuous thrust for first body
    if i == 1
        rv_dot(4:6,i) = rv_dot(4:6,i) + controls(rv(:,i),t);
    end
end

end

function uvect = rk4(mu,u0,dt,T)
% 4-Stage Runge-Kutta Propagator

uvect = repmat(u0,1,1,T/dt);

for jj = 1 : T/dt
    
    uvect(:,:,jj) = u0;
    tvect = linspace(dt*(jj-1),dt*jj,4);
    
    % ============= RK4 =============
    k1 = dt * f( u0,mu,tvect(1));
    k2 = dt * f((u0 + k1/2),mu,tvect(2));
    k3 = dt * f((u0 + k2/2),mu,tvect(3));
    k4 = dt * f((u0 + k3)  ,mu,tvect(4));
    u0 = u0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end

end