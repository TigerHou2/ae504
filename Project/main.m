close all
clear;clc

mu = 1.327e20;
r0 = [150.63e9, 0, 0 ;...
      217.10e9, 0, 0 ]';
v0 = [0, 29.72e3, 0 ;...
      0, 24.13e3, 0]';
x0 = [r0;v0];

dt = 100;
tf = 2500000;

mu = 0;

% run simulation
xvect = dynamics(mu,x0,dt,tf);

% assign trajectory colors
cmap = hsv(size(xvect,2));

% plot trajectory
figure(1)
hold on
scatter3(0,0,0,36,'Filled','MarkerFaceColor','k')
for i = 1:size(xvect,2)
    x = xvect(1,i,:);
    y = xvect(2,i,:);
    z = xvect(3,i,:);
    plot3(x(:),y(:),z(:),'MarkerSize',16,'Color',cmap(i,:))
end

% plot target
scatter3(-217.10e9,0,0,36,'Filled','MarkerFaceColor','r')

hold off

axis equal
% 
% % plot thrust vector
% figure(2)
% tvect = linspace(0,tf,tf/dt+1);
% thrust_vect = nan(3,size(uvect,3));
% for i = 1:tf/dt
%     thrust_vect(:,i) = controls(uvect(:,1,i),tvect(i));
% end
% plot3(thrust_vect(1,:),thrust_vect(2,:),thrust_vect(3,:))