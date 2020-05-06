close all
clear;clc

r0 = [0 0]';
v0 = [0 0]';
rf = [10 10]';
vf = [1 -1]';

umax = 2;

g = [0 -1]';

dt = 0.01;

t0 = 0;
tf = 10;

%% Controls Calc
syms c1 c2 n1 n2 tf L0 real

c  = [c1 c2]';
lr = [n1 n2]';
lv = @(t) lr*t + c;

lv0 = lv(0);
lvf = lv(tf);

% L0 = 1;

u0 = -umax * lv0 / sqrt(lv0'*lv0);
uf = -umax * lvf / sqrt(lvf'*lvf);

eqn1 = 0 == L0 + lr'*v0 + lv0'*(g+u0);
eqn2 = 0 == L0 + lr'*vf + lvf'*(g+uf);
eqn3 = vf-v0 == int(lv,0,tf);

[c1,c2,n1,n2,tf,L0] = vpasolve([eqn1,eqn2,eqn3],[c1,c2,n1,n2,tf,L0]);
c1 = double(c1);
c2 = double(c2);
n1 = double(n1);
n2 = double(n2);
tf = double(tf);

%% Simulation

tf = 5;

tt = 0:dt:tf;
rr = nan(2,length(tt));
vv = nan(2,length(tt));
uu = nan(2,length(tt));

r = r0;
v = v0;

for i = 1:length(tt)
%     u = [n1 n2]'*tt(i) + [c1 c2]';
    u = [-0.2 -0.69]'*tt(i) + [1.09 3.02]';
    u = u / norm(u) * umax;
    r = r + dt * v;
    v = v + dt * (g + u);
    rr(:,i) = r;
    vv(:,i) = v;
    uu(:,i) = u;
end

figure(1)
hold on
plot(rr(1,:),rr(2,:))
scatter(rf(1),rf(2))
xlabel('x')
ylabel('y')
title('Trajectory')
figure(2)
plot(uu(1,:),uu(2,:))
xlabel('x')
ylabel('y')
title('Controls')