function xvect = dynamics_v2(mu,x0,dt,t0,tf,lr,lv)
%DYNAMICS_V2 Summary of this function goes here
%   Detailed explanation goes here

xvect = repmat(x0,1,1,(tf-t0)/dt);

pp = [-lv;lr];

for jj = 1 : (tf-t0)/dt
    
    pu = pp(4:6);
    
    disp(norm(pu))
    
    if norm(pu)>1
        pu = pu / norm(pu) * 0.02;
    else
        pu = zeros(3,1);
    end
    
    xvect(:,:,jj) = x0;
    
    % ============= RK4 =============
    k1 = dt * f( x0,mu,pu);
    k2 = dt * f((x0 + k1/2),mu,pu);
    k3 = dt * f((x0 + k2/2),mu,pu);
    k4 = dt * f((x0 + k3)  ,mu,pu);
    x0 = x0 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    % ============= Primer Propagation ==============
    r = x0(1:3);
    G = mu/norm(r)^5*(3*(r*r')-(r'*r)*eye(3));
    pp = expm([[zeros(3), eye(3)];[G,zeros(3)]]*(dt))*pp;
    
end

end


%% function definitions
function [rv_dot,u] = f(rv,mu,u)
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
        rv_dot(4:6,i) = rv_dot(4:6,i) + u;
    end
    % ----- end of control function -----
end

end