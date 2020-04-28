close all
clear;clc

mu = 1;
r0 = [1 0 0]';
v0 = [0 1 0]';
u0 = [r0;v0];

dt = 0.001;
tf = 50;

% run simulation
uvect = dynamics(mu,u0,dt,tf);

% assign trajectory colors
cmap = hsv(size(uvect,2));

% plot trajectory
figure(1)
hold on
scatter3(0,0,0,36,'Filled','MarkerFaceColor','k')
for i = 1:size(uvect,2)
    x = uvect(1,i,:);
    y = uvect(2,i,:);
    z = uvect(3,i,:);
    plot3(x(:),y(:),z(:),'MarkerSize',10,'Color',cmap(i,:))
end
hold off

axis equal

% plot thrust vector
figure(2)
tvect = linspace(0,tf,tf/dt+1);
thrust_vect = nan(3,size(uvect,3));
for i = 1:tf/dt
    thrust_vect(:,i) = controls(uvect(:,1,i),tvect(i));
end
plot3(thrust_vect(1,:),thrust_vect(2,:),thrust_vect(3,:))