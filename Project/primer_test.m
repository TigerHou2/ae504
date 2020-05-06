close all
clear;clc

syms p1 p2 p3 k real

pa = [p1 p2 p3]';

mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 10]';

dt = 1e-3;

umax = 0.3;

r1 = r0 + dt * v0;
v1 = v0 + dt * (-mu * r0 / norm(r0)^3 + pa/sqrt(pa'*pa)*umax);

pb = pa + dt * (-v0 - 3 * dt * (-mu * r0 / norm(r0)^3));

eqn1 = k*pa == 2*r0;
eqn2 = pa'*pa == 1;
eqn3 = k*pb == 2*r1 - 3*dt*v1;

soln = solve([eqn1,eqn2,eqn3],[p1,p2,p3,k])