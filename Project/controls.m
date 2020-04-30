function u = controls(rv,mu,uvect,tvect,t)
%CONTROLS Summary of this function goes here
%   Detailed explanation goes here

%% V4
umax = 1;
u = nan(3,1);
u(1) = sign(tvect(1)-t)*uvect(1);
u(2) = sign(tvect(2)-t)*uvect(2);
u(3) = sign(tvect(3)-t)*uvect(3);

u = u + mu * rv(1:3) / norm(rv(1:3))^3;

if norm(u) > umax
    disp('Warning: thrust exceeds tolerance!')
end

if t > tvect(4)
    u = zeros(3,1);
    return
end

% if (mod(t,50000)<1)
%     hold on
%     quiver3(rv(1),rv(2),rv(3),u(1)*8e10,u(2)*8e10,u(3)*8e10,'b');
%     hold off
% end

%% V3
% tgt = [-217.10e9, 0, 0, 0, -24.13e3, 0]';
% umax = 0.3;
% 
% r1 = tgt(1)-rv(1);
% r2 = tgt(2)-rv(2);
% r3 = tgt(3)-rv(3);
% 
% Lv1 = uprev(1);
% Lv2 = uprev(2);
% Lv3 = uprev(3);
% 
% u = -[...
%     (mu*(2*Lv1*r1^2 + 3*Lv2*r1*r2 + 3*Lv3*r1*r3 - Lv1*r2^2 - Lv1*r3^2))...
%     /(r1^2 + r2^2 + r3^2)^(5/2);...
%     (mu*(2*Lv1*r1^2 + 3*Lv2*r1*r2 + 3*Lv3*r1*r3 - Lv1*r2^2 - Lv1*r3^2))...
%     /(r1^2 + r2^2 + r3^2)^(5/2);...
%     (mu*(2*Lv1*r1^2 + 3*Lv2*r1*r2 + 3*Lv3*r1*r3 - Lv1*r2^2 - Lv1*r3^2))...
%     /(r1^2 + r2^2 + r3^2)^(5/2)];
% 
% u = u / norm(u) * umax;
% 
% if (mod(t,50000)<1)
%     hold on
%     quiver3(rv(1),rv(2),rv(3),u(1)*8e10,u(2)*8e10,u(3)*8e10,'b');
%     hold off
% end

%% V2
% tgt = [-217.10e9, 0, 0, 0, -24.13e3, 0]';
% umax = 0.1;
% 
% rw = 0.7e-28;
% vw = 1;
% dpos = tgt(1:3)-rv(1:3);
% dvel = tgt(4:6)-rv(4:6);
% 
% u = rw*dpos*norm(dpos) + vw*dvel/(norm(dpos)+1);
% 
% if (norm(dpos)<=1e5) && (norm(dvel)<=10)
%     u = zeros(3,1);
%     disp("Hol' up")
%     return
% else
%     u = u / norm(u) * umax;
% end
% 
% if (mod(t,100000)<1)
%     hold on
%     quiver3(rv(1),rv(2),rv(3),u(1)*5e10,u(2)*5e10,u(3)*5e10,'b');
%     hold off
% end

%% v1
% u = -0.05 * rv(1:3) / norm(rv(1:3)) + 0.05 * sin(t/3);

end

