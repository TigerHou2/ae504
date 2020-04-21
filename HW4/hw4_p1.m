%%

close all
clear;clc

close all
clear;clc

% initial position
r0 = [-10;50];
v0 = [0;0];

% landing sites
tar1 = [100;0];
tar2 = [150;0];

% control bounds
bounds = [-2,2];

% controls

vx_m1 = 2;  % 2
vy_m1 = 4.5;  % 3

vx_m2 = 3.5;  % 3
vy_m2 =-3.3;  % -1

state_ini = [r0;v0];
state_m1  = [10;terrain(10);vx_m1;vy_m1] + [-0.4;1;0;0];
state_m2  = [50;terrain(50);vx_m2;vy_m2] + [ 1;0;0;0];
state_end = [tar1;0;0];

% control points

[u1,t1] = controls_v2(   0   ,bounds,state_ini,state_m1);
[u2,t2] = controls_v2(max(t1),bounds,state_m1 ,state_m2);
[u3,t3] = controls_v2(max(t2),bounds,state_m2 ,state_end);

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
hold off

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

%%
close all
clear;clc

% initial position
r0 = [-10;50];
v0 = [0;0];

% landing sites
tar1 = [100;0];
tar2 = [150;0];

% control bounds
bounds = [-2 2;
          -2 2];

% controls
state_ini = [r0;v0];
state_end = [tar1;0;0];

% control points
t = linspace(0,52,100000);

u = controls(t,bounds,state_ini,state_end);

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
hold off

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
min(r(2,:))
disp(['Control cost:' num2str(trapz(t,abs(u(1,:)))+trapz(t,abs(u(2,:))))])