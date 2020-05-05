close all
clear;clc

t = linspace(0,1,10000);
dt = t(2) - t(1);

ee = exp(1);

x = 0;
cost = 0;

for i = 1:length(t)-1
%     u = 1/(ee-1)*exp(t(i));
    u = 1;
%     u = 2*(1-t(i));
    x = x + u*dt;
    cost = cost + dt * u^2;
end

disp(['  xf = ' num2str(x)])
disp(['cost = ' num2str(cost)])

%%
close all
clear;clc

syms C t real

nu = -2*C-4;

u = -0.5*nu*t - 0.5*C;

int(u^2,t,0,1)