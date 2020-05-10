close all
clear;clc

r0 = [0 0]';
v0 = [0 1]';
rf = [6 6]';
vf = [-1 -1]';
umax = 2;
g = [0 -1]';
dt = 0.01;
t0 = 0;

params = {r0,v0,rf,vf,umax,g,dt,t0};

%% Controls Calc

X0 = [-0.07,-1.07,0.002,0.05,10];
options = optimoptions('fsolve','Display','iter'...
                               ,'PlotFcn','optimplotx');
xx = fsolve(@(X) trajSolver(X,params),X0,options);

c1 = xx(1);
c2 = xx(2);
n1 = xx(3);
n2 = xx(4);
tf = xx(5);

%% Simulation

tt = 0:dt:tf;
rr = nan(2,length(tt));
vv = nan(2,length(tt));
uu = nan(2,length(tt));

r = r0;
v = v0;

for i = 1:length(tt)
    u = [n1 n2]'*tt(i) + [c1 c2]';
    % getting rid of u=0 singularity
    if norm(u) < 1e-6
        u = [0;0];
    else
        u = - u / norm(u) * umax;
    end
    % getting rid of fsolve precision error
    u(abs(u)<1e-3) = 0;
    r = r + dt * v;
    v = v + dt * (g + u);
    rr(:,i) = r;
    vv(:,i) = v;
    uu(:,i) = u;
end

figure(2)
hold on
plot(rr(1,:),rr(2,:))
scatter(rf(1),rf(2))
xlabel('x')
ylabel('y')
title('Trajectory')
figure(3)
plot(uu(1,:),uu(2,:))
xlabel('x')
ylabel('y')
title('Controls')

%% function definitions

function F = trajSolver(X,params)

    r0 = params{1};
    v0 = params{2};
    rf = params{3};
    vf = params{4};
    umax = params{5};
    g = params{6};

    c1 = X(1);
    c2 = X(2);
    n1 = X(3);
    n2 = X(4);
    tf = X(5);
    
    L0 = 0.1;

    c  = [c1 c2]';
    lr = [n1 n2]';
    lv = @(t) (lr*t + c);
    u  = @(t) -umax * (lr*t + c) / sqrt((lr*t + c)'*(lr*t + c));
    v  = @(t) v0 + integral(u,0,t,'ArrayValued',true) + g*t;

    F(1) =  L0 + lr'*v0 + lv(0)'*(u(0)+g);
    F(2:3) =  -vf + v0 + integral(u,0,tf,'ArrayValued',true);
    F(4:5) =  -rf + r0 + integral(v,0,tf,'ArrayValued',true);

end