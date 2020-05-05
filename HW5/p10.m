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
Ldot = [-diff(H,x1);...
        -diff(H,x2);...
        -diff(H,x3)];

% H = 0
simplify(subs(H,u,uu)) == 0

%%
close all
clear;clc

A = [1 2 7; ...
     2 2 4; ...
     1 3 3];
b = [1 3 3]';

Z = b*b';

x = [1 0 0]';

% --------------
% syms u x1 x2 x3 real
% xx = [x1 x2 x3]';
% eqn = xx' * ( A' * (2*Z*xx - 2*Z*(A*xx+b*u)) + xx ) == u^2;
% pretty(isolate(eqn,u))
% --------------

t = linspace(0,10,10000);
dt = t(2) - t(1);

Kk = lqr(A,b,eye(3),1,zeros(3,1));
Q = eye(3);
R = 1;
P = riccati(A,b,Q,R);
K = inv(R)*(b'*P);

cost = 0;

for i = 1:length(t)
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
%     u = (35921*x1^2 + 122154*x1*x2 + 201156*x1*x3 + 103854*x2^2 ...
%        + 342038*x2*x3 + 281625*x3^2)^(1/2) - 323*x2 - 532*x3 - 190*x1;
%     u = -dot(b,x);
    u = -K*x;
    x = x + dt * ( A*x + b*u );
    cost = cost + dt * (x'*x + u^2);
end

%% function definitions
function P = riccati(A,B,Q,R)
    P = Q;
    while true
        Pp = P + 0.001 * (A'*P + P*A - P*B*inv(R)*B'*P + Q);
        if max(abs(Pp-P))<1e-10
            P = Pp;
            return
        else
            P = Pp;
        end
    end
end