%% Simulation Initialization

close all
clear;clc

% define spacecraft and physics
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 10]';
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path

% define simulation parameters
dt = 1000;
t0 = 0;
tf = 10000;

% define target's current state and controls
tgt = [-217.10e9, 0, 0, ... % position
       0, -24.13e3, 10]';   % velocity
t_ini = 0;
t_end = 10000;
x0 = [x0,tgt];

%% Simulation Execution

% calculate target state at end time
[re,ve] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);

% create symbolic vector of Lagrange multipliers
syms Lr1 Lr2 Lr3 Lv1 Lv2 Lv3 real
Lr = [Lr1;Lr2;Lr3];
Lv = [Lv1;Lv2;Lv3];

xvect = dynamics_v2(mu,x0,dt,t0,tf,Lr,Lv);

% assign trajectory colors
cmap = hsv(size(xvect,2));
