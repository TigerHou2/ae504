function a = dynamics(u)
%DYNAMICS Summary of this function goes here
%   Detailed explanation goes here

g = 1.62;

a = [u(1,:);u(2,:)-g];

end

