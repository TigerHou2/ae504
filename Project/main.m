%% Simulation Initialization

close all
clear;clc

% define spacecraft and physics
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 10]';
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path

% define simulation parameters
dt = 2000;
t0 = 0;
tf = 50000000;

% define target's current state and controls
tgt = [217.10e9, 0, 100e9, ... % position
       0, 24.13e3, 0]';   % velocity
t_ini = 0;
t_end = 15000000;
x0 = [x0,tgt];

% define number of thrust segments
segments = 1;

%% Simulation Execution

% calculate target state at end time
[re,ve] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);

% split time into segments
tvect = linspace(t_ini,t_end,segments+1)';

% interpolate orbital elements between start target at segment times
% --- first, get orbital parameters
[a_ini,e_ini,i_ini,omg_ini,w_ini,f_ini] = Get_Orb_Params(r0(:,1),v0(:,1),mu);
[a_end,e_end,i_end,omg_end,w_end,f_end] = Get_Orb_Params(re(:,1),ve(:,1),mu);
% --- next, interpolate orbital parameters (currently we will assume same
% orbit direction for origin and target.
params_ini = [a_ini,e_ini,i_ini,omg_ini,w_ini,f_ini];
params_end = [a_end,e_end,i_end,omg_end,w_end,f_end];
params_vect = repmat(params_ini,segments+1,1) ...
            + linspace(0,1,segments+1)'*(params_end-params_ini);
% --- create end points for flight segments
xx = nan(6,segments+1);
xx(:,1) = x0(:,1);
xx(:,end) = [re;ve];
% --- create empty output vectors
x_curr = x0;
xvect = [];
dv = [];

% run simulation
for i = 1:segments+1
    % --- account for transition from 360 -> 0
    pt = params_vect(i,:);
    pt = double(pt<0).*[0,0,0,2*pi,2*pi,2*pi] + pt;
    params_vect(i,:) = pt;
    % --- compute interpolated r,v vectors
    if i ~= 1 && i ~= segments+1
        [r,v] = Get_Orb_Vects(pt,mu);
        xx(:,i) = [r;v];
        hold on
        scatter3(r(1),r(2),r(3),32,'magenta','Filled')
        hold off
    end
    % --- compute control for current segment
    if i ~= 1
        [uc,tc] = init_control_v1(xx(:,i-1),xx(:,i),tvect(i-1),tvect(i));
        if i == segments+1
            [xvect_curr,dv_curr] = dynamics(mu,x_curr,dt,tvect(i-1),tf,uc,[tc;tvect(i)]);
        else
            [xvect_curr,dv_curr] = dynamics(mu,x_curr,dt,tvect(i-1),tvect(i),uc,[tc;tvect(i)]);
        end
        % --- update bodies in sim
        x_curr = xvect_curr(:,:,end);
        xvect = cat(3,xvect,xvect_curr);
        dv = [dv; dv_curr];
    end
end

% assign trajectory colors
cmap = hsv(size(xvect,2));

%% Rendering

% plot trajectory
figure(1)
hold on
scatter3(0,0,0,32,'Filled','MarkerFaceColor','k')
for i = 1:size(xvect,2)
    x = xvect(1,i,:);
    y = xvect(2,i,:);
    z = xvect(3,i,:);
    plot3(x(:),y(:),z(:),'LineWidth',1.2,'Color',cmap(i,:))
end

% plot target
scatter3(re(1),re(2),re(3),32,cmap(end,:),'Filled')
% plot all bodies
scatter3(xvect(1,:,end),xvect(2,:,end),xvect(3,:,end),32,cmap,'Filled')

hold off
axis equal
view([-1,-1,1])

%% Output

idx = (t_end-t0)/dt;

% display position and velocity
disp(['Position: ' num2str(xvect(1:3,1,idx)') ' m'])
disp(['Velocity: ' num2str(xvect(4:6,1,idx)') ' m/s'])
disp([' Delta-v: ' num2str(sum(dv)) ' m/s'])