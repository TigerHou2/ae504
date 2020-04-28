function accel = controls(rv,t)
%CONTROLS Summary of this function goes here
%   Detailed explanation goes here

accel = -0.05 * rv(1:3) / norm(rv(1:3)) + 0.05 * sin(t/3);

end

