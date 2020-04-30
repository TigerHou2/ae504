close all
clear;clc

% define all bodies (except target)
mu = 1.327e20;
r0 = [150.63e9, 0, 0 ;... % s/c
      150.63e9, 0, 0]';
v0 = [0, 29.72e3, 0 ;... % s/c
      0, 29.72e3, 0]';
x0 = [r0;v0];

% define simulation parameters
dt = 1000;
t0 = 0;
tf = 20000000;

% define target's current state and controls
tgt = [217.10e9, 0, 0, ... % position
       0, 24.13e3, 0]';   % velocity
t_ini = 0;
t_end = 10000000;
[re,ve] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24); % tgt future state
[u,t] = init_control_v1(x0(:,1),[re;ve],t_ini,t_end);
t(end+1) = t_end;

x0 = [x0,tgt];

% run simulation
xvect = dynamics(mu,x0,dt,t0,tf,u,t);

% assign trajectory colors
cmap = hsv(size(xvect,2));

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

idx = (t_end-t0)/dt;

% display position and velocity
disp(['Position: ' num2str(xvect(1:3,1,idx)')])
disp(['Velocity: ' num2str(xvect(4:6,1,idx)')])