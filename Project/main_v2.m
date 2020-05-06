%% Simulation Initialization

close all
clear;clc

% define spacecraft and physics
mu = 1.327e20;
r0 = [150.63e9, 0, 0]';
v0 = [0, 29.72e3, 0]';
x0 = [[r0,r0];[v0,v0]]; % duplicate to show original path

% define simulation parameters
dt = 1000;
t0 = 0;
tf = 21000000;

% define target's current state and controls
tgt = [-217.10e9, 0, -100e9, ... % position
       0, -24.13e3, 0]';   % velocity
t_ini = 0;
t_end = 20000000;
x0 = [x0,tgt];

% define number of thrust segments
segments = 1;

%% Simulation Execution

% calculate target state at end time
[rf,vf] = TimeProp_V3(tgt(1:3),tgt(4:6),mu,(t_end-t_ini)/3600/24);

% split time and weights into segments
tvect = linspace(t_ini,t_end,segments+1)';
weights = linspace(0,1,segments+1);
xx = x0(:,1);

% interpolate h vectors
hi = cross(r0(:,1),v0(:,1));
hf = cross(rf,vf);
hvect = [ linspace(hi(1),hf(1),segments+1) ;...
          linspace(hi(2),hf(2),segments+1) ;...
          linspace(hi(3),hf(3),segments+1) ];

% interpolate r magnitudes
ri_norm = norm(r0(:,1));
rf_norm = norm(rf(:,1));
rnorm_vect = linspace(ri_norm,rf_norm,segments+1);

% find r vectors
% --- for each h, project r0 and rf onto the plane defined by h
r0_proj = r0(:,1) - (hvect ./ vecnorm(hvect,2,1)) .* ((r0(:,1)'*hvect) ./ vecnorm(hvect,2,1));
rf_proj = rf(:,1) - (hvect ./ vecnorm(hvect,2,1)) .* ((rf(:,1)'*hvect) ./ vecnorm(hvect,2,1));
% --- rotate r0 in the direction of h toward rf by an angle
%     defined by the number of segments for interpolation
angles = linspace(0,1,segments+1)...
       .*acos(sum(r0_proj.*rf_proj)./vecnorm(r0_proj,2,1)...
                                   ./vecnorm(rf_proj,2,1));
rr = nan(size(r0_proj));
for i = 1:size(rr,2)
    rr(:,i) = rotv(r0_proj(:,i),hvect(:,i),angles(i));
end
% --- this vector is then scaled between r0 and rf by a factor
%     defined by the number of segments
rr = rr ./ vecnorm(rr,2,1) .* rnorm_vect;
% --- finally, the velocity is defined by the position vector and h
vv = nan(size(rr));
for i = 2:size(vv,2)
    vv(:,i) = cross(hvect(:,i),rr(:,i)) / (rr(:,i)'*rr(:,i));
%     vv(:,i) = norm(vv(:,i)) * (rr(:,i)-rr(:,i-1)) / norm(rr(:,i)-rr(:,i-1));
end

% create empty output vectors
x_curr = x0;
xvect = [];
dv = [];

%% run simulation
for i = 1:segments+1
    % --- compute interpolated r,v vectors
    if i ~= segments+1
        rv = [rr(:,i+1);vv(:,i+1)];
        hold on
        scatter3(rv(1),rv(2),rv(3),32,'magenta','Filled')
        hold off
    end
    % --- compute control for current segment
    if i ~= 1
        [uc,tc] = init_control_v1(xx,rv,tvect(i-1),tvect(i));
        if i == segments+1
            [xvect_curr,dv_curr] = dynamics(mu,x_curr,dt,tvect(i-1),tf,uc,[tc;tvect(i)]);
        else
            [xvect_curr,dv_curr] = dynamics(mu,x_curr,dt,tvect(i-1),tvect(i),uc,[tc;tvect(i)]);
        end
        % --- update bodies in sim
        x_curr = xvect_curr(:,:,end);
        xx = x_curr(:,1);
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
scatter3(rf(1),rf(2),rf(3),32,cmap(end,:),'Filled')
% plot all bodies
scatter3(xvect(1,:,end),xvect(2,:,end),xvect(3,:,end),32,cmap,'Filled')

hold off
axis equal
view([-1 -1 1])

%% Output

idx = (t_end-t0)/dt;

% display position and velocity
disp(['Position: ' num2str(xvect(1:3,1,idx)') ' m'])
disp(['Velocity: ' num2str(xvect(4:6,1,idx)') ' m/s'])
disp([' Delta-v: ' num2str(sum(dv)) ' m/s'])