close all
clear;clc

syms x1 x2 x3 L1 L2 L3 u real

A = [1 2 7; ...
     2 2 4; ...
     1 3 3];
b = [1 3 3]';

x = [x1 x2 x3]';
lambda = [L1 L2 L3]';

% Hamiltonian
H = x'*x + u^2 + lambda'*A*x + lambda'*b*u;

% dH/du = 0
uu = rhs( isolate(diff(H,u)==0,u) );

% dH/dx = Ldot
Ldot = [diff(H,x1);...
        diff(H,x2);...
        diff(H,x3)];

% H = 0
simplify(subs(H,u,uu)) == 0