%% gravity_correction.m
%
% Author: 
%   Tiger Hou
%
% Description:
%   This script calculates the controls for a transfer between 
%   two planets by assuming there are no body forces and then
%   compensating for gravity in the controls.
%

%% Simulation Initialization

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
tgt = [217.10e9, 0    , 100e9, ... % position
       0    , 24.13e3 , 0    ]';   % velocity
t_ini = t0;
t_end = 10e6;

    % % case 2
    % tgt = [217.10e9, 0    , 0, ... % position
    %        0    , 24.13e3 , 0    ]';   % velocity
    % t_ini = t0;
    % t_end = 25.7e6;

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
params = {r0,v0,rf,vf,umax,scale};

%% Solve Optimal Control

% solve for optimal control policy (without gravity)
X0 = [0.108,-0.889,-0.283,2.414e-7,-2.044e-6,-6.526e-7,8.857e5] ./ scale;
options = optimoptions('fsolve','Display','iter'...
                               ,'PlotFcn','optimplotx');
xx = fsolve(@(X) trajSolver(X,params),X0,options);
xx = xx .* scale;
c1 = xx(1);
c2 = xx(2);
c3 = xx(3);
n1 = xx(4);
n2 = xx(5);
n3 = xx(6);
tf = xx(7);
p = @(t) ( -[n1 n2 n3]'*t + [c1 c2 c3]' ) * (t<=tf); % lv

%% Run Simulation

% run simulation
[xvect,uvect,tvect] = run_sim(mu,x0,p,t_final,dt,umax);

% assign trajectory colors
cmap = hsv(size(xvect,2));

%% Rendering

% plot trajectory
figure(2)
hold on
scatter3(0,0,0,32,'Filled','MarkerFaceColor','k')
pp = gobjects(1,size(xvect,2));
for i = 1:size(xvect,2)
    x = xvect(1,i,:);
    y = xvect(2,i,:);
    z = xvect(3,i,:);
    pp(i) = plot3(x(:),y(:),z(:),'LineWidth',1.2,'Color',cmap(i,:));
end

% plot target
pt = scatter3(rf(1),rf(2),rf(3),32,cmap(end,:));
% plot all bodies
% scatter3(xvect(1,:,end),xvect(2,:,end),xvect(3,:,end),32,cmap,'Filled')

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

% plot controls
figure(3)
hold on
tt = tvect;
uu = uvect;
tt(:,vecnorm(uu,2,1)==0) = [];
uu(:,vecnorm(uu,2,1)==0) = [];
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

% plot controls by axis
figure(4)
hold on
plot(tt,uu(1,:),'LineWidth',1.2)
plot(tt,uu(2,:),'LineWidth',1.2)
plot(tt,uu(3,:),'LineWidth',1.2)
hold off
xlabel('time, $s$');
ylabel('$\Gamma u ,\; m/s^2$')
grid(gca,'minor')
grid on
legend('$\Gamma u_x$','$\Gamma u_y$','$\Gamma u_z$')
latexify(19,13)

%% Output

% display position and velocity
disp(['Position: ' num2str(xvect(1:3,1,end)') ' m'])
disp(['Velocity: ' num2str(xvect(4:6,1,end)') ' m/s'])
disp([' Delta-v: ' num2str(trapz(linspace(0,tf,size(uvect,2)),...
                                 vecnorm(uvect,2,1))) ' m/s'])

%% Functions -- Optimal Control

function F = trajSolver(X,params)

    r0 = params{1};
    v0 = params{2};
    rf = params{3};
    vf = params{4};
    umax = params{5};
    scale = params{6};

    xx = X .* scale;
    c1 = xx(1);
    c2 = xx(2);
    c3 = xx(3);
    n1 = xx(4);
    n2 = xx(5);
    n3 = xx(6);
    tf = xx(7);
    
    L0 = 1;

    c  = [c1 c2 c3]';
    lr = [n1 n2 n3]';
    lv = @(t) (-lr*t + c);
    u  = @(t) -umax * (-lr*t + c) / sqrt((-lr*t + c)'*(-lr*t + c));
    v  = @(t) v0 + integral(u,0,t,'ArrayValued',true);

    F(1) =  L0 + lr'*v0 + lv(0)'*(u(0));
    F(2:4) =  -vf + v0 + integral(u,0,tf,'ArrayValued',true);
    F(5:7) =  -rf + r0 + integral(v,0,tf,'ArrayValued',true);
    
    F(2:4) = F(2:4) * norm(rf) / norm(vf);

end

%% Functions -- Simulation

function [xvect,uvect,tvect] = run_sim(mu,x0,p,t_final,dt,umax)

    steps = ceil(t_final/dt);
    dt = t_final/steps;
    xvect = nan(6,size(x0,2),steps-1);
    tvect = linspace(0,t_final,steps);
    uvect = nan(3,steps-1);
    rv = x0;

    for jj = 1 : steps-1
        
        % remove u = 0 singularity
        if norm(p(tvect(jj))) < 1e-6
            u = zeros(3,1);
        else
            u = -p(tvect(jj)) / norm(p(tvect(jj))) * umax;
        end
        % remove fsolve precision error
        u(abs(u)<1e-3) = 0;
        
        % reapply gravity
        if norm(p(tvect(jj))) == 0
            % control ends, do not apply gravity correction
        else
            u = u + mu * rv(1:3,1) / norm(rv(1:3,1))^3;
        end
        % store uvect
        uvect(:,jj) = u;

        xvect(:,:,jj) = rv;

        % ============= RK4 =============
        k1 = dt * f( rv,mu,u);
        k2 = dt * f((rv + k1/2),mu,u);
        k3 = dt * f((rv + k2/2),mu,u);
        k4 = dt * f((rv + k3)  ,mu,u);
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    end
    
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