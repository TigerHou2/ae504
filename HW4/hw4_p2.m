%% 

close all
clear;clc

syms A B C real

tarx = 100;

eqn1 = 50 == 100*A - 10*B + C;
eqn2 = 0  == tarx^2*A + tarx*B + C;
eqn3 = 400/3 + 1 == 2500*A + 50*B + C;

soln = solve([eqn1,eqn2,eqn3],[A,B,C]);
a = double(soln.A);
b = double(soln.B);
c = double(soln.C);

x = linspace(-10,tarx,10000);
para = @(xx) a.*(xx.^2) + b.*xx + c;
y = para(x);

peak = [-b/2/a,para(-b/2/a)];

plot(x,y,'LineWidth',1)
hold on
terrain_x = linspace(-20,160,5000);
terrain_y = terrain(terrain_x);
plot(terrain_x,terrain_y,'LineWidth',1);
hold off

ux = x(2:end)-x(1:end-1);
uy = y(2:end)-y(1:end-1);

%%
close all
clc

% initial position
r0 = [-10;50];
v0 = [0;0];

% landing sites
tar1 = [100;0];
tar2 = [150;0];

% control bounds
bounds = [-2,2];

% controls

g = 1.62;

% left eject point on parabola
rx1 = 20;
ry1 = para(rx1);
vx1 = (peak(1)-rx1)/sqrt(2*(peak(2)-ry1)/g);
vy1 = sqrt(2*(peak(2)-ry1)*g);

% right eject point on parabola
rx2 = 60;
ry2 = para(rx2);
vx2 = (rx2-peak(1))/sqrt(2*(peak(2)-ry2)/g);
vy2 = -sqrt(2*(peak(2)-ry2)*g);

state_ini = [r0;v0];
state_m1  = [rx1;ry1;vx1;vy1];
state_m2  = [rx2;ry2;vx2;vy2];
state_end = [tar1;0;0];

% control points

[u1,t1] = controls_v3(   0   ,bounds,state_ini,state_m1);
t2 = linspace(t1(end),t1(end)+(rx2-rx1)/vx1,20000);
u2 = zeros(2,length(t2));
[u3,t3] = controls_v3(max(t2),bounds,state_m2 ,state_end);

u1(:,end) = [];
t1(end) = []; % get rid of redundant point overlapping t1 and t2
u2(:,end) = [];
t2(end) = [];

u = [u1,u2,u3];
t = [t1,t2,t3];

% responses
a = dynamics(u);

v = v0 + cumtrapz(t,a,2);
r = r0 + cumtrapz(t,v,2);

% plot terrain
terrain_x = linspace(-20,160,5000);
terrain_y = terrain(terrain_x);
plot(terrain_x,terrain_y,'LineWidth',1);

% plot trajectory
hold on
plot(r(1,:),r(2,:),'LineWidth',1);
plot(x,y,'LineWidth',1)
scatter([rx1,rx2],[ry1,ry2],'Filled')
hold off

legend('Terrain','Lander Trajectory','Parabolic Trajectory','Injection Points')

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','ticklabelinterpreter'),...
                            'ticklabelinterpreter','latex')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set( gca, 'Color', [1 1 1] )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [30 12])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0.2 1.2 30 12])
set(gcf, 'PaperPosition', [0.2 1.2 30 12])
grid(gca,'minor')
grid on

% output target offset
disp(['Distance from first target: ' newline ...
      '     x = ' num2str(r(1,end)-tar1(1)) newline ...
      '     y = ' num2str(r(2,end)-tar1(2))])
disp(['Distance from second target: ' newline ...
      '     x = ' num2str(r(1,end)-tar2(1)) newline ...
      '     y = ' num2str(r(2,end)-tar2(2))])
disp(['Final velocity:' newline ...
      '    vx = ' num2str(v(1,end)) newline ...
      '    vy = ' num2str(v(2,end))])
disp(['Control cost:' num2str(trapz(t,abs(u(1,:)))+trapz(t,abs(u(2,:))))])
disp(['Flight time: ' num2str(max(t))])