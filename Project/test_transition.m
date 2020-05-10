%% setup
% ----------------------------
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660016262.pdf
%

close all
clear;clc

% define spacecraft and physics
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 0]';

% define simulation parameters
dt = 300;
t0 = 0;
t_final = 35e6;

% define target state
tgt = [217.10e9, 0, 100e9, ... % position
       0, 24.13e3, 0]';   % velocity
t_ini = t0;
t_end = 10e6;

% define control limits
umax = 1;

% set up objects for plotting
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path
x0 = [x0,tgt];

% calculate target state at end time
[rf,vf] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);
xf = [rf;vf];

% make life easier for fsolve by scaling variables
scale = [1e1,1e1,1e0,1e-6,1e-6,1e-7,1e6];

% gather parameters for solving controls
params = {r0,v0,rf,vf,mu,umax,scale};

%% solve for controls
% ----------------------------

% gravity gradient matrix
G0 = mu/norm(r0)^5*(3*r0*r0'-norm(r0)^2*eye(3));

% transition matrix
phi = @(dt) [ [ eye(3) + G0*dt^2/2, eye(3)*dt + G0*dt^3/6 ];...
              [ G0*dt             , eye(3) + G0*dt^2/2] ];

% initial guess
% X0 = [-0.108,0.889,0.283,2.414e-7,-2.044e-6,-6.526e-7,8.8e5] ./ scale;
X0 = [-0.6131,0.7008,0.2128,1.2118e-6,-1.3903e-6,-4.2032e-7,1.0350e6] ./ scale;
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',700);
xx = fsolve(@(X) trajSolver(X,params,phi), X0, options);

xx  = xx .* scale;
p0  = xx(1:3)';
pd0 = xx(4:6)';
tf  = xx(7);

pp0 = [p0;pd0];

p = @(t) [eye(3) zeros(3)] * (phi(t) * pp0);
u = @(t) (norm(p(t))>1e-3) * umax * p(t) / sqrt(p(t)'*p(t)) * (t<tf);


%% execute simulation
% ----------------------------

% create initial state vector
rv0 = x0;

% simulate control and gather position data
[rv_vect,uvect] = run_sim(mu,rv0,u,dt,t_final);

% assign trajectory colors
cmap = hsv(size(rv_vect,2));

%% plotting

% plot trajectory
figure(2)
hold on
scatter3(0,0,0,32,'Filled','MarkerFaceColor','k')
for i = 1:size(rv_vect,2)
    x = rv_vect(1,i,:);
    y = rv_vect(2,i,:);
    z = rv_vect(3,i,:);
    plot3(x(:),y(:),z(:),'LineWidth',1.2,'Color',cmap(i,:))
end

% plot target
scatter3(rf(1),rf(2),rf(3),32,cmap(end,:))
% plot all bodies
scatter3(rv_vect(1,:,end),rv_vect(2,:,end),rv_vect(3,:,end),32,cmap,'Filled')

hold off
axis equal
view([-1,-1,1])
xlabel('x')
ylabel('y')
zlabel('z')
grid(gca,'minor')
grid on
latexify(19,13)

% plot controls
uu = uvect;
uu(:,vecnorm(uu,2,1)==0) = [];
figure(3)
hold on
plot3(uu(1,:),uu(2,:),uu(3,:),'LineWidth',1.2)
scatter3(uu(1,end),uu(2,end),uu(3,end))
hold off
axis equal
view([-1,-1,1])
xlabel('x')
ylabel('y')
zlabel('z')
grid(gca,'minor')
grid on
latexify(19,13)

%% trajectory solving functions

function F = trajSolver(X,params,phi)

    r0 = params{1};
    v0 = params{2};
    rf = params{3};
    vf = params{4};
    mu = params{5};
    umax = params{6};
    scale = params{7};
    
    xx = X .* scale;
    n1  = xx(1);
    n2  = xx(2);
    n3  = xx(3);
    nd1 = xx(4);
    nd2 = xx(5);
    nd3 = xx(6);
    tf  = xx(7);
    
    p0  = [n1 n2 n3]';      % -lv
    pd0 = [nd1 nd2 nd3]';   %  lr
    pp0 = [p0;pd0];
    
    p = @(t) [eye(3) zeros(3)] * (phi(t) * pp0);
    u = @(t) (norm(p(t))>1e-3) * umax * p(t) / sqrt(p(t)'*p(t));
    
    F(1) = 1 + pd0'*v0 - p0'*(-mu*r0/norm(r0)^3+u(0));
    F(2:7) = -[rf;vf] + [r0;v0] + propagate(r0,v0,mu,u,tf);
    
    F(5:7) = F(5:7) * norm(rf) / norm(vf);
    
end

function [rvf] = propagate(r0,v0,mu,u,tf)
    
    tvect = linspace(0,tf,501);
    dt = tvect(2) - tvect(1);
    rv = [r0;v0];
    for i = 1:length(tvect)-1
        
        k1 = dt * f( rv,mu,u(tvect(i)));
        k2 = dt * f((rv + k1/2),mu,u(tvect(i)));
        k3 = dt * f((rv + k2/2),mu,u(tvect(i)));
        k4 = dt * f((rv + k3)  ,mu,u(tvect(i)));
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
        
    end
    rvf = rv;

end

%% simulation functions

function [rv_vect,uvect,tvect] = run_sim(mu,rv0,u,delta_t,t_final)

    steps = ceil(t_final/delta_t);
    tvect = linspace(0,t_final,steps+1);
    dt = tvect(2) - tvect(1);
    
    rv_vect = nan(6,size(rv0,2),steps);
    uvect  = nan(3,steps);
    
    rv = rv0;

    for jj = 1 : steps
    
        uvect(:,jj) = u(tvect(jj));

        rv_vect(:,:,jj) = rv;

        % ============= RK4 =============
        k1 = dt * f( rv,mu,u(tvect(jj)));
        k2 = dt * f((rv + k1/2),mu,u(tvect(jj)));
        k3 = dt * f((rv + k2/2),mu,u(tvect(jj)));
        k4 = dt * f((rv + k3)  ,mu,u(tvect(jj)));
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    end
    
    % trim tvect
    tvect(end) = [];

end

function [rv_dot,u] = f(rv,mu,u)
% Two-Body Gravity

    N = size(rv,2);
    rv_dot = zeros(size(rv));

    for i = 1:N
        rv_dot(1:3,i) = rv(4:6,i); % vel
        rv_dot(4:6,i) = - mu * rv(1:3,i) / norm(rv(1:3,i))^3; % accel

        if i == 1
            rv_dot(4:6,i) = rv_dot(4:6,i) + u; % control
        end
    end

end