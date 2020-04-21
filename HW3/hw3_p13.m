close all
clear;clc

x0 = [1;0];

A = [ 2 0;
      1 2];
B = [1;0];

Q = [1 0;
     0 1];
Qn = Q;
R = 1;

d = [3;3];

N = 3; % number of states from initial to end

P = zeros(size(Q,1),size(Q,2),N+1);
P(:,:,N+1) = Qn;

F = zeros(size(B,2),size(B,1),N+1);

for i = 1:N
    
    P1 = P(:,:,N-i+2);
    P(:,:,N-i+1) = A' * P1 * A...
               - A' * P1 * B * inv(R+B'*P1*B)*B'*P1*A + Q;
    F(:,:,N-i+1) = -inv(R+B'*P1*B)*B'*P1*A;
    
end

T = N;
% [reward,u,xn] = V(x0,[],A,B,Q,R,T,N,d)

coeff = 0.578;
[reward,u,xn] = V(x0,[],A,B,Q,R,F,T,N,d,coeff)

%% function definitions

function [reward,u,xn] = V(x,u,A,B,Q,R,F,T,N,d,coeff)
    u(end+1) = F(:,:,N+1-T)*(x+coeff*d);
    uu = u(end);
    xn = A*x + B*uu + d;
    if T == 0
        reward = x'*Q*x;
        u(end) = [];
    else
        [reward,u] = V(xn,u,A,B,Q,R,F,T-1,N,d,coeff);
        reward = reward + x'*Q*x + uu'*R*uu;
    end
end

% function [reward,u,xn] = V(x,u,A,B,Q,R,T,N,d)
%     u(end+1) = -2*inv(2*R+1*B'*Q*B)*B'*Q*(A*x+d);
%     uu = u(end);
%     xn = A*x + B*uu + d;
%     if T == 0
%         reward = x'*Q*x;
%         u(end) = [];
%     else
%         [reward,u,xn] = V(xn,u,A,B,Q,R,T-1,N,d);
%         reward = reward + x'*Q*x + uu'*R*uu;
%     end
% end