%% optimal_control.m
%
% Author: 
%   Tiger Hou
%
% Description:
%   This script optimizes the controls for a continuous high-thrust
%   spacecraft performing a brachistochrone transfer between two planets.
%

%% setup
% ----------------------------

close all
clear;clc

% define spacecraft origin
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 0]';

% define RK4 simulation parameters
dt = 200;
t0 = 0;
t_final = 35e6;

% define control limit
umax = 1;

% define target and calculate target state at desired time
tgt = [217.10e9, 0    , 100e9, ... % position
       0    , 24.13e3 , 0    ]';   % velocity
t_ini = t0;
t_end = 10e6;
[rf,vf] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);
xf = [rf;vf];

    % % case 2
    % tgt = [217.10e9, 0    , 0, ... % position
    %        0    , 24.13e3 , 0    ]';   % velocity
    % t_ini = t0;
    % t_end = 25.7e6;

% set up objects for plotting
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path of spacecraft
x0 = [x0,tgt];          % append target to list of objects

% improve fsolve by scaling variables
% to be roughly the same order of magnitude
scale = [1e1,1e1,1e0,1e-6,1e-6,1e-7,1e6];

% gather parameters for solving optimal control problem
params = {r0,v0,rf,vf,mu,umax,scale};

%% solve for controls
% ----------------------------

% gravity gradient matrix
G0 = mu/norm(r0)^5*(3*r0*r0'-norm(r0)^2*eye(3));

% truncated transition matrix
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660016262.pdf
phi = @(dt) [ [ eye(3) + G0*dt^2/2, eye(3)*dt + G0*dt^3/6 ];...
              [ G0*dt             , eye(3) + G0*dt^2/2] ];

% initial guess
% modify initial guess by scaling
X0 = [-0.1084,0.8880,0.2842,2.4495e-7,-2.0423e-6,-6.5465e-7,8.8275e5] ./ scale;
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',700);
xx = fsolve(@(X) trajSolver(X,params,phi), X0, options);

% restore output by undoing scale factor
xx  = xx .* scale;
p0  = xx(1:3)';
pd0 = xx(4:6)';
tf  = xx(7);

pp0 = [p0;pd0];

% define primer vector p(t) and control input u(t)
% both functions output 3-by-1 vectors
p = @(t) [eye(3) zeros(3)] * (phi(t) * pp0);
u = @(t) (norm(p(t))>1e-5) * umax * p(t) / sqrt(p(t)'*p(t)) * (t<tf);


%% execute simulation
% ----------------------------

% simulate control and gather position data
[rv_vect,uvect,tvect] = run_sim(mu,x0,u,dt,tf,t_final);

% assign trajectory colors
cmap = hsv(size(rv_vect,2));

%% plotting
% ----------------------------

% plot trajectory
figure(2)
hold on
scatter3(0,0,0,32,'Filled','MarkerFaceColor','k')
pp = gobjects(1,size(rv_vect,2));
for i = 1:size(rv_vect,2)
    x = rv_vect(1,i,:);
    y = rv_vect(2,i,:);
    z = rv_vect(3,i,:);
    pp(i) = plot3(x(:),y(:),z(:),'LineWidth',1.2,'Color',cmap(i,:));
end

% plot target
pt = scatter3(rf(1),rf(2),rf(3),32,cmap(end,:));

hold off
axis equal
view([-1,-1,1])
xlabel('x, $m$')
ylabel('y, $m$')
zlabel('z, $m$')
grid(gca,'minor')
grid on
legend([pp,pt],'Spacecraft','Origin','Destination','Target','Location','Best')
latexify(19,13)

% plot controls, ignore u = 0 after tf
tt = tvect;
uu = uvect;
tt(:,vecnorm(uu,2,1)==0) = [];
uu(:,vecnorm(uu,2,1)==0) = [];
figure(3)
hold on
plot3(uu(1,:),uu(2,:),uu(3,:),'LineWidth',1.2)
scatter3(uu(1,end),uu(2,end),uu(3,end))
% --- select points along control curve to show path of u
idx = 2175;
indices = [idx-8;idx-6;idx-4;idx-2;idx;idx+2;idx+4;idx+6;idx+8];
num_arrows = length(indices);
for i = 1:num_arrows
    Ov(1,i) = uu(1,indices(i));
    Ov(2,i) = uu(2,indices(i));
    Ov(3,i) = uu(3,indices(i));
    Dv(1,i) = uu(1,indices(i)+1) - uu(1,indices(i)-1);
    Dv(2,i) = uu(2,indices(i)+1) - uu(2,indices(i)-1);
    Dv(3,i) = uu(3,indices(i)+1) - uu(3,indices(i)-1);
end
Dv = Dv ./ vecnorm(Dv,2,1);
quiver3(Ov(1,:),Ov(2,:),Ov(3,:),Dv(1,:),Dv(2,:),Dv(3,:),...
        'MaxHeadSize',1,'LineWidth',1,'AutoScaleFactor',0.7)
hold off
axis equal
view([-1,-1,1])
xlabel('$\Gamma u_x ,\; m/s^2$')
ylabel('$\Gamma u_y ,\; m/s^2$')
zlabel('$\Gamma u_z ,\; m/s^2$')
grid(gca,'minor')
grid on
latexify(19,15)

% plot controls by axis
figure(4)
hold on
plot(tt,uu(1,:),'LineWidth',1.2)
plot(tt,uu(2,:),'LineWidth',1.2)
plot(tt,uu(3,:),'LineWidth',1.2)
hold off
xlabel('time, $s$');
ylabel('$\Gamma u ,\; m/s^2$')
legend('$\Gamma u_x$','$\Gamma u_y$','$\Gamma u_z$')
grid(gca,'minor')
grid on
latexify(19,13)

%% trajectory solving functions
% ----------------------------------
% these functions are used to numerically integrate the trajectory
% and then solve for the optimal control

function F = trajSolver(X,params,phi)
% Boundary Conditions for Optimal Trajectory

    % collect B.C.s and I.C.s of the problem formulation
    r0 = params{1};
    v0 = params{2};
    rf = params{3};
    vf = params{4};
    mu = params{5};
    umax = params{6};
    scale = params{7};
    
    % collect optimization variables
    xx = X .* scale; % apply scaling factor to make fsolve work better
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
    
    % define primer vector p(t) and control input u(t)
    % both functions output 3-by-1 vectors
    p = @(t) [eye(3) zeros(3)] * (phi(t) * pp0);
    u = @(t) (norm(p(t))>1e-5) * umax * p(t) / sqrt(p(t)'*p(t));
    
    % boundary conditions aresatisfied when F = 0
    F(1) = 1 + pd0'*v0 - p0'*(-mu*r0/norm(r0)^3+u(0));
    F(2:7) = -[rf;vf] + propagate(r0,v0,mu,u,tf);
    
    % apply scaling factor to F
    % so velocity and position are weighed equally
    F(5:7) = F(5:7) * norm(rf) / norm(vf);
    
end

function [rvf] = propagate(r0,v0,mu,u,tf)
% Simulate Brachistochrone Trajectory
%   r0 : 3-by-1 position vector of spacecraft
%   v0 : 3-by-1 velocoty vector of spacecraft
%   mu : gravitational parameter of central body
%   u  : anonymous control function u(t) with time argument
    
    tvect = linspace(0,tf,3001);
    dt = tvect(2) - tvect(1);
    rv = [r0;v0];
    for i = 1:length(tvect)-1
        
        % ============= RK4 =============
        k1 = dt * f( rv,mu,u(tvect(i)));
        k2 = dt * f((rv + k1/2),mu,u(tvect(i)));
        k3 = dt * f((rv + k2/2),mu,u(tvect(i)));
        k4 = dt * f((rv + k3)  ,mu,u(tvect(i)));
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
        
    end
    rvf = rv;

end

%% simulation functions
% ----------------------------------
% these functions are used to simulate the orbit
% after the optimal control has been calculated

function [rv_vect,uvect,tvect] = run_sim(mu,rv0,u,delta_t,tf,t_final)
% Simulate Brachistochrone Trajectory
%   mu     : gravitational parameter of central body
%   rv0    : 6-by-N vector of initial position, velocity of N bodies
%            we assume that the first column is the spacecraft
%   u      : anonymous control function u(t) with time argument
%   delta_t: RK4 time step
%   tf     : control input final time
%   t_final: simulation final time

    % propagate through time using two step sizes:
    % --- while thrusting, use a small time step delta_t
    %     this is because the gravity vector changes quickly
    steps_ctrl  = ceil(tf/delta_t);
    tvect_ctrl  = linspace(0,tf,steps_ctrl+1);
    tvect_ctrl(end)  = []; % trim end of tvect_ctrl
    % --- while coasting, use a larger time step delta_t * 10
    %     because the gravity vector now changes slowly
    steps_coast = ceil((t_final-tf)/delta_t/10);
    tvect_coast = linspace(tf,t_final,steps_coast+1);
    % --- finally, combine time vectors
    tvect = [tvect_ctrl, tvect_coast];
    
    steps = length(tvect)-1;
    
    % allocate memory to store state vectors and controls
    rv_vect = nan(6,size(rv0,2),steps);
    uvect  = nan(3,steps);
    
    % we will use the variable rv to step through RK4
    rv = rv0;

    for jj = 1 : steps
    
        % record controls and state
        uvect(:,jj) = u(tvect(jj));
        rv_vect(:,:,jj) = rv;
        
        % calculate dt at the current step
        dt = tvect(jj+1) - tvect(jj);

        % ============= RK4 =============
        k1 = dt * f( rv        ,mu,u(tvect(jj)));
        k2 = dt * f((rv + k1/2),mu,u(tvect(jj)));
        k3 = dt * f((rv + k2/2),mu,u(tvect(jj)));
        k4 = dt * f((rv + k3)  ,mu,u(tvect(jj)));
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    end
    
    % trim tvect because the final time is not used in the calculations
    tvect(end) = [];

end

function [rv_dot,u] = f(rv,mu,u)
% Restricted Two-Body Gravity Model
%   rv : 6-by-N vector of position, velocity of N bodies
%        we assume that the first column is the spacecraft
%   mu : gravitational parameter of central body
%   u  : 3-by-1 control input (acceleration)

    N = size(rv,2);
    rv_dot = zeros(size(rv));

    % iterate through N orbiting bodies
    for i = 1:N
        rv_dot(1:3,i) = rv(4:6,i); % vel
        rv_dot(4:6,i) = - mu * rv(1:3,i) / norm(rv(1:3,i))^3; % accel

        if i == 1
            rv_dot(4:6,i) = rv_dot(4:6,i) + u; % control
        end
    end

end