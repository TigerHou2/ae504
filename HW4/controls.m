function u = controls(t,bounds,state_ini,state_end)
%CONTROLS Summary of this function goes here
%   Detailed explanation goes here

u = nan(2,length(t));

tmax = max(t);
xmin = min(bounds(1,:));
xmax = max(bounds(1,:));
ymin = min(bounds(2,:));
ymax = max(bounds(2,:));

% case-specific inputs
h_start = 0.50558;
h_reverse = 1-(1-h_start)/2;

% horizontal control
    xnodes = [15.1674,22.5837,30];
    xctrls = [0,xmax,xmin,0];

    % -- start and end conditions
    xcond{length(xctrls)} = t>=xnodes(end);
    xcond{1} = t<xnodes(1);

    % -- internal conditions
    for i = 2:length(xnodes)
        xcond{i} = t>=xnodes(i-1) & t<xnodes(i);
    end

    % -- apply conditions to control
    for i = 1:length(xctrls)
        u(1,xcond{i}) = xctrls(i);
    end

% vertical control
    ynodes = [21,25.94];
    yctrls = [ymax,ymin,ymax];

    % -- start and end conditions
    ycond{length(yctrls)} = t>=ynodes(end);
    ycond{1} = t<ynodes(1);

    % -- internal conditions
    for i = 2:length(ynodes)
        ycond{i} = t>=ynodes(i-1) & t<ynodes(i);
    end

    % -- apply conditions to control
    for i = 1:length(yctrls)
        u(2,ycond{i}) = yctrls(i);
    end

end

