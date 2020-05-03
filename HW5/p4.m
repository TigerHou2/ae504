close all
clear;clc

syms x(t)

ode = diff(x,t,2) == 2*x^3+x;

soln(t) = dsolve(ode)