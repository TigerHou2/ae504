close all
clear;clc

syms ls li S I alpha beta real
syms I positive

H = 10*I + alpha + beta + ls*(-(1-alpha)*I*S) ...
                        + li*( (1-alpha)*I*S - beta*I);

eqn1 = diff(H,alpha) == 0;
eqn2 = diff(H,beta)  == 0;

soln = solve([eqn1,eqn2],[ls,li]);

Ls = soln.ls;
Li = soln.li;

lsd = simplify(-diff(H,S));
lid = simplify(-diff(H,I));

Lsd = simplify(subs(lsd,[ls,li],[Ls,Li]));
Lid = simplify(subs(lid,[ls,li],[Ls,Li]));

disp(['Ls     =  ' char(Ls)])
disp(['Li     =  ' char(Li)])
disp(['Ls_dot =  ' char(Lsd)])
disp(['Li_dot =  ' char(Lid)])

sd = -(1-alpha)*I*S;
id =  (1-alpha)*I*S - beta*I;

eqn_ls = (I*S*sd - (S-1)*(id*S+sd*I)) / (I^2*S^2) == id;
eqn_li = -10*I^2 + I*(alpha+beta-1) == id;

eqn_alpha_1 = isolate(simplify(eqn_ls),alpha);
eqn_alpha_2 = isolate(simplify(eqn_li),alpha);

[a,b] = solve([eqn_alpha_1,eqn_alpha_2],[alpha,beta]);

%% simulation

close all
clear;clc

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

ls = (ss-1)/ii/ss;
li = 1/ii;

cost = 0;

for k = 1:length(t)

    signA =  sign(ls)*sign(li)*sign(li/ii+ls/ss);
    signB = -sign(li);

    A = max(signA,0);
    B = max(signB,0);

    ls = ls + dt * ( (A-1)/ss );
    li = li + dt * ( (A+B-1)/ii - 10 );

%     ls = ls + dt * ( ii*(li-ls)*(A-1) );
%     li = li + dt * ( li*(B+ss*(A-1))-ss*ls*(A-1)-10 );
    
    avect(k) = A;
    bvect(k) = B;
    
    ss = ss - dt * ( (1-A)*ii*ss );
    ii = ii + dt * ( (1-A)*ii*ss - B*ii );
    rr = rr + dt * ( B*ii );
    
    svect(k) = ss;
    ivect(k) = ii;
    rvect(k) = rr;
    
    cost = cost + ( 10*ii + A + B ) * dt;
    
end

figure(1)

subplot(2,1,1)
hold on
plot(t,avect,'--r','LineWidth',1.2)
plot(t,bvect,'--b','LineWidth',1.2)
hold off
legend('a','b','Location','Best')
grid on

subplot(2,1,2)
hold on
% plot(t,svect)
plot(t,ivect,'LineWidth',1.2)
plot(t,rvect,'LineWidth',1.2)
hold off
legend('Infected','Recovered','Location','Best')
grid on

sgtitle(['Cost = ' num2str(cost)])
% sgtitle('No Control Demonstration')
latexify

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [20 10])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0.2 1.2 20 10])
set(gcf, 'PaperPosition', [0.2 1.2 20 10])