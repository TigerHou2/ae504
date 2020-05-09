%% setup
% ----------------------------
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660016262.pdf
%

close all
clear;clc

% define spacecraft and physics
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 10]';
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path

% define simulation parameters
dt = 2000;
t0 = 0;

% define target state
tgt = [0, 217.10e9, 20e9, ... % position
      -24.13e3, 0, 0]';   % velocity
  
% define control limits
umax = 1;

% define fsolve variable scale factors
scale = [1e1,1e1,1e0,1e-6,1e-6,1e-7,1e6];

% pass parameters to trajectory solver
params = {r0,v0,tgt(1:3),tgt(4:6),mu,umax,scale};

%% solve for controls
% ----------------------------

G0 = mu/norm(r0)^5*(3*r0*r0'-norm(r0)^2*eye(3));
phi = @(dt) [ [ eye(3) + G0*dt^2/2, eye(3)*dt + G0*dt^3/6 ];...
              [ G0*dt             , eye(3) + G0*dt^2/2] ];

X0 = [0.642,-0.8259,-0.0822,-1.271e-6,-1.658e-6,1.64e-7,1.01e6] ./ scale;
options = optimoptions('fsolve','Display','iter'...
                               ,'PlotFcn','optimplotx'...
                               ,'MaxFunctionEvaluations',500);
xx = fsolve(@(X) trajSolver(X,params,phi), X0, options);

xx = xx .* scale;
lr = xx(4:6)';
lv = xx(1:3)';
tf = xx(7);


%% execute simulation
% ----------------------------

% create initial state vector
rv0 = [[r0;v0],tgt];
% simulate and gather position data
[rv_vect,primer] = run_sim(mu,rv0,lr,lv,tf,dt,umax);
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
scatter3(tgt(1),tgt(2),tgt(3),32,cmap(end,:),'Filled')
% plot all bodies
scatter3(rv_vect(1,:,end),rv_vect(2,:,end),rv_vect(3,:,end),32,cmap,'Filled')

hold off
axis equal
view([-1,-1,1])

% plot primer norm
figure(3)
tvect = linspace(0,tf,size(primer,2));
plot(tvect,vecnorm(primer,2,1))

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
    
    p0  = [n1 n2 n3]';
    pd0 = [nd1 nd2 nd3]';
    pp0  = [p0;pd0];
    
    p = @(t) [eye(3) zeros(3)] * (phi(t) * pp0);
    
    u = @(t) (norm(p(t))>1) * umax * p(t) / sqrt(p(t)'*p(t)) + zeros(3,1);
    
    F(1) = -1 + p0'*p0;
    F(2:7) = -[rf;vf] + [r0;v0] + propagate(r0,v0,mu,u,tf);
    
    F(5:7) = F(5:7) * norm(rf) / norm(vf);
    
end

function [rvf] = propagate(r0,v0,mu,u,tf)

    tvect = linspace(0,tf,2001);
    dt = tvect(2) - tvect(1);
    r = r0;
    v = v0;
    for i = 1:length(tvect)-1
        rtemp = r + v * dt;
        v = v + u(tvect(i)) - mu*r/norm(r)^3;
        r = rtemp;
    end
    rf = r;
    vf = v;
    rvf = [rf;vf];

end

%% simulation functions

function [rv_vect,primer] = run_sim(mu,rv0,lr,lv,tf,dt,umax)

    steps = ceil(tf/dt);
    dt = tf/steps;
    rv_vect = nan(6,size(rv0,2),steps-1);
    primer  = nan(3,steps-1);
    pp = [-lv;lr];
    rv = rv0;

    for jj = 1 : steps-1
    
        pu = pp(4:6);
        primer(:,jj) = pu;
        if norm(pu)>1
            pu = pu / (pu'*pu) * umax;
        else
            pu = pu * 0;
        end

        rv_vect(:,:,jj) = rv;

        % ============= RK4 =============
        k1 = dt * f( rv,mu,pu);
        k2 = dt * f((rv + k1/2),mu,pu);
        k3 = dt * f((rv + k2/2),mu,pu);
        k4 = dt * f((rv + k3)  ,mu,pu);
        rv = rv + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

        % ============= Primer Propagation ==============
        r = rv(1:3,1);
        G = mu/(r'*r)^5*(3*(r*r')-(r'*r)*eye(3));
        pp = pp + [[zeros(3), eye(3)];[G,zeros(3)]]*pp*dt;

    end

end

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