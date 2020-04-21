close all
clear;clc

x0 = [1;0];
A = [2 0; 1 2];
B = [1;0];
d = [3;3];

syms u0 u1 u2 real

u = [u0;u1;u2];

x1 = A*x0 + B*u(1) + d;
x2 = A*x1 + B*u(2) + d;
x3 = A*x2 + B*u(3) + d;

objfun = x0'*x0 + x1'*x1 + x2'*x2 + x3'*x3 + u'*u;

eqn1 = diff(objfun,u0) == 0;
eqn2 = diff(objfun,u1) == 0;
eqn3 = diff(objfun,u2) == 0;

soln = solve([eqn1,eqn2,eqn3],u);
U0 = double(soln.u0);
U1 = double(soln.u1);
U2 = double(soln.u2);

cost = double(subs(objfun,u',[U0,U1,U2]))