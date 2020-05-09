%%
close all
clear;clc

syms r1 r2 r3 v1 v2 v3 Lr1 Lr2 Lr3 Lv1 Lv2 Lv3 real
syms mu real
syms L(t)

Lr = [Lr1 Lr2 Lr3]';
Lv = [Lv1 Lv2 Lv3]';
r = [r1 r2 r3]';
v = [v1 v2 v3]';

lambda = [Lr;Lv];
x = [r;v];

H = Lr'*v + Lv'*(mu*r/norm(r)^3);

Ldot = -[diff(H,r1);...
         diff(H,r2);...
         diff(H,r3);...
         diff(H,v1);...
         diff(H,v2);...
         diff(H,v3)];
simplify(Ldot(1))

%% no gravity testing
close all
clear;clc

r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 0]';
% rt = [-217.10e9, 0, 0]';
% vt = [0, -24.13e3, 0]';
rt = [0, 217.10e9, 0]';
vt = [-24.13e3, 0, 0]';
umax = 1;
tf = 2000000;

% r0 = [0,0,0]';
% v0 = [0,0,0]';
% rt = [4,0,0]';
% vt = [0,0,0]';
% tf = 4;

syms x y z real
syms tx ty tz positive
eqn1 = v0(1) + x*tx - x*(tf-tx) - vt(1) == 0;
eqn2 = v0(2) + y*ty - y*(tf-ty) - vt(2) == 0;
eqn3 = v0(3) + z*tz - z*(tf-tz) - vt(3) == 0;
eqn4 = r0(1) + v0(1)*tx + 0.5*x*tx^2 + (v0(1)+x*tx)*(tf-tx) - 0.5*x*(tf-tx)^2 - rt(1) == 0;
eqn5 = r0(2) + v0(2)*ty + 0.5*y*ty^2 + (v0(2)+y*ty)*(tf-ty) - 0.5*y*(tf-ty)^2 - rt(2) == 0;
eqn6 = r0(3) + v0(3)*tz + 0.5*z*tz^2 + (v0(3)+z*tz)*(tf-tz) - 0.5*z*(tf-tz)^2 - rt(3) == 0;
% eqn4 = r0(1) + v0(1)*tf + 0.5*x*(2*tx-tf) - x*tx*(tx-tf) - rt(1) == 0;
% eqn5 = r0(2) + v0(2)*tf + 0.5*y*(2*ty-tf) - y*ty*(ty-tf) - rt(2) == 0;
% eqn6 = r0(3) + v0(3)*tf + 0.5*z*(2*tz-tf) - z*tz*(tz-tf) - rt(3) == 0;

soln1 = solve([eqn1,eqn4],[x,tx]);
soln2 = solve([eqn2,eqn5],[y,ty]);
soln3 = solve([eqn3,eqn6],[z,tz]);

u = nan(3,1);
t = nan(3,1);
u(1) = double(soln1.x(1));
u(2) = double(soln2.y(1));
u(3) = double(soln3.z(1));
t(1) = double(soln1.tx(1));
t(2) = double(soln2.ty(1));
t(3) = double(soln3.tz(1));

u
t

x = u(1);
y = u(2);
z = u(3);
tx = t(1);
ty = t(2);
tz = t(3);