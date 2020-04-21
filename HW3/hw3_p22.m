close all
clear;clc

x0 = [0; 10; 0; 0];

A = [ 0.994  0.026  0     -32.2;
     -0.094  0.376  820    0;
      0     -0.002  0.332  0;
      0      0      1      1];
B = [0; -32.7; -2.08; 0];

Q = [1 0 0 0;
     0 1 0 0;
     0 0 0 0;
     0 0 0 1];
Qn = Q;
R = 1;

N = 10;

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
[reward,u,xn] = V(x0,[],A,B,Q,R,F,T,N)

%% function definitions

function [reward,u,xn] = V(x,u,A,B,Q,R,F,T,N)
    u(end+1) = F(:,:,N+1-T)*x;
    uu = u(end);
    xn = A*x + B*uu;
    if T == 0
        reward = -x'*Q*x;
        u(end) = [];
    else
        [reward,u] = V(xn,u,A,B,Q,R,F,T-1,N);
        reward = reward - x'*Q*x - uu'*R*uu;
    end
end