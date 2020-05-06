%% Simulation Initialization

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
tf = 15000000;

% define target's current state and controls
tgt = [217.10e9, 0, 100e9, ... % position
       0, 24.13e3, 0]';   % velocity
t_ini = 0;
t_end = 15000000;
x0 = [x0,tgt];

%% Simulation Execution

% calculate target state at end time
[re,ve] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);
xf = [re;ve];

% run simulation
[xvect,dv] = dynamics_v3(mu,x0,dt,t0,tf);

% assign trajectory colors
cmap = hsv(size(xvect,2));

%% Rendering

% plot trajectory
figure(1)
hold on
scatter3(0,0,0,32,'Filled','MarkerFaceColor','k')
for i = 1:size(xvect,2)
    x = xvect(1,i,:);
    y = xvect(2,i,:);
    z = xvect(3,i,:);
    plot3(x(:),y(:),z(:),'LineWidth',1.2,'Color',cmap(i,:))
end

% plot target
scatter3(re(1),re(2),re(3),32,cmap(end,:),'Filled')
% plot all bodies
scatter3(xvect(1,:,end),xvect(2,:,end),xvect(3,:,end),32,cmap,'Filled')

hold off
axis equal
view([-1,-1,1])

%% Output

idx = (t_end-t0)/dt;

% display position and velocity
disp(['Position: ' num2str(xvect(1:3,1,idx)') ' m'])
disp(['Velocity: ' num2str(xvect(4:6,1,idx)') ' m/s'])
disp([' Delta-v: ' num2str(sum(dv)) ' m/s'])