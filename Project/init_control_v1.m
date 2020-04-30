function [u,t] = init_control_v1(rv0,rvf,t0,tf)
%INIT_CONTROL_V1 Summary of this function goes here
%   Detailed explanation goes here

r0 = rv0(1:3);
v0 = rv0(4:6);
rt = rvf(1:3);
vt = rvf(4:6);

syms x y z real
syms tx ty tz positive
eqn1 = v0(1) + x*(tx-t0) - x*(tf-tx) - vt(1) == 0;
eqn2 = v0(2) + y*(ty-t0) - y*(tf-ty) - vt(2) == 0;
eqn3 = v0(3) + z*(tz-t0) - z*(tf-tz) - vt(3) == 0;
eqn4 = r0(1) + v0(1)*(tx-t0) + 0.5*x*(tx-t0)^2 + (v0(1)+x*(tx-t0))*(tf-tx) - 0.5*x*(tf-tx)^2 - rt(1) == 0;
eqn5 = r0(2) + v0(2)*(ty-t0) + 0.5*y*(ty-t0)^2 + (v0(2)+y*(ty-t0))*(tf-ty) - 0.5*y*(tf-ty)^2 - rt(2) == 0;
eqn6 = r0(3) + v0(3)*(tz-t0) + 0.5*z*(tz-t0)^2 + (v0(3)+z*(tz-t0))*(tf-tz) - 0.5*z*(tf-tz)^2 - rt(3) == 0;

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

end

