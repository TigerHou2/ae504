close all
clear;clc

t = linspace(0,1,100000);
dt = t(2) - t(1);

x = 0;
cost = 0;
A = 0;

xx = nan(length(t)-1,1);

for i = 1:length(t)-1
%     u = sign(0.5-t(i))*0.4;
    nu = -1.08531866;
    C = -nu/2;
    u = (nu*t(i)+C)*sqrt(1/(1-(nu*t(i)+C)^2));
%     u = sign(0.6105-t(i)) * 0.25*exp(t(i));
    x = x + u * dt;
    cost = cost + dt * sqrt(1+u^2);
    A = A + x*dt;
    xx(i) = x;
end

disp(['  xf = ' num2str(x)])
disp(['cost = ' num2str(cost)])
disp(['   A = ' num2str(A)])

plot(t(1:end-1),xx)
axis equal

%%
close all
clear;clc

syms nu t tt real

C = -nu/2;

u = (nu*t+C)*sqrt(1/(1-(nu*t+C)^2));

eqn = int(int(u,t,0,tt),tt,0,1) == 0.1;
Nu = vpasolve(eqn,nu)