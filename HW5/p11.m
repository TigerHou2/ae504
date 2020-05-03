close all
clear;clc

syms ls li S I alpha beta real
syms I positive

eqn1 = 1 + (ls-li)*I*S == 0;
eqn2 = 1 - li*I == 0;

soln = solve([eqn1,eqn2],[ls,li]);

ls = soln.ls;
li = soln.li;

lsd = simplify(-(li-ls)*(1-alpha)*I);
lid = simplify(-(li-ls)*(1-alpha)*S + li*beta - 10);

disp(['Ls     =  ' char(ls)])
disp(['Li     =  ' char(li)])
disp(['Ls_dot =  ' char(lsd)])
disp(['Li_dot =  ' char(lid)])

sd = -(1-alpha)*I*S;
id =  (1-alpha)*I*S - beta*I;

eqn_ls = (I*S*sd - (S-1)*(id*S+sd*I)) / (I^2*S^2) == id;
eqn_li = -10*I^2 + I*(alpha+beta-1) == id;

eqn_alpha_1 = isolate(simplify(eqn_ls),alpha);
eqn_alpha_2 = isolate(simplify(eqn_li),alpha);

[a,b] = solve([eqn_alpha_1,eqn_alpha_2],[alpha,beta]);

%% simulation
t = linspace(0,10,1000);
dt = t(2) - t(1);

ss = 0.999;
ii = 0.001;
rr = 0;

avect = nan(size(t));
bvect = nan(size(t));

svect = nan(size(t));
ivect = nan(size(t));
rvect = nan(size(t));

cost = 0;

for k = 1:length(t)
    
%     A = double(subs(a,[S,I],[ss,ii]));
%     B = double(subs(b,[S,I],[ss,ii]));

    A = 0.27;
    B = 0.27;
    
    A = max(-1,min(1,A));
    B = max(-1,min(1,B));
    
    avect(k) = A;
    bvect(k) = B;
    
    ss = ss - ( (1-A)*ii*ss ) * dt;
    ii = ii + ( (1-A)*ii*ss - B*ii ) * dt;
    rr = rr + ( B*ii ) * dt;
    
    svect(k) = ss;
    ivect(k) = ii;
    rvect(k) = rr;
    
    cost = cost + ( 10*ii + A + B ) * dt;
    
end

% figure(1)
% hold on
% plot(t,avect,'--r','LineWidth',1)
% plot(t,bvect,'--b','LineWidth',1)
% hold off
% legend('a','b','Location','Best')
% title(['Cost = ' num2str(cost)])

figure(2)
hold on
% plot(t,svect)
plot(t,ivect)
plot(t,rvect)
hold off
legend('Infected','Recovered','Location','Best')
title(['Cost = ' num2str(cost)])