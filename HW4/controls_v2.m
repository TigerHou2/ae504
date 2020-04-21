function [u,t] = controls_v2(tMin,bounds,state_ini,state_end)
%CONTROLS_V2 Summary of this function goes here
%   Detailed explanation goes here

res = 10000;
tmin = 0;

g = 1.62;

umin = min(bounds);
umax = max(bounds);

% set up symbolic variables for single switch
syms tt ttx tty tmax ratio positive
    
% vertical motion

    if state_end(2) >= state_ini(2)
        % vehicle ascends
        yctrls = [umax,0];
    else
        % vehicle descends
        yctrls = [0,umax];
    end

    u1 = yctrls(1);
    v1(tt) = state_ini(4) + int(u1-g,tt,tmin,tt);
    r1(tt) = state_ini(2) + int(v1,tt,tmin,tt);
    vp = v1(tty);
    rp = r1(tty);
    
    u2 = yctrls(2);
    v2(tt) = vp + int(u2-g,tt,tt,tmax);
    r2(tt) = rp + int(v2,tt,tt,tmax);
    vf = v2(tty);
    rf = r2(tty);
    
    eqn_y_1 = rf == state_end(2);
    eqn_y_2 = vf == state_end(4);
    eqn_y_3 = tmin <= tty;
    eqn_y_4 = tty <= tmax;
    
% solve for time required

soln = solve([eqn_y_1,eqn_y_2,eqn_y_3,eqn_y_4],[tty,tmax]);
ty = double(soln.tty);
tMax = double(soln.tmax);
        
% horizontal motion

    if state_end(1) >= state_ini(1)
        % right translation
        xctrls = [umax,umin];
    else
        % left translation
        xctrls = [umin,umax];
    end
    xCtrls = xctrls * ratio;
    
    u1 = xCtrls(1);
    v1(tt) = state_ini(3) + int(u1,tt,tmin,tt);
    r1(tt) = state_ini(1) + int(v1,tt,tmin,tt);
    vp = v1(ttx);
    rp = r1(ttx);
    
    u2 = xCtrls(2);
    v2(tt) = vp + int(u2,tt,tt,tMax);
    r2(tt) = rp + int(v2,tt,tt,tMax);
    vf = v2(ttx);
    rf = r2(ttx);
    
    eqn_x_1 = rf == state_end(1);
    eqn_x_2 = vf == state_end(3);
    eqn_x_3 = tmin <= ttx;
    eqn_x_4 = ttx <= tMax;
        
soln = solve([eqn_x_1,eqn_x_2,eqn_x_3,eqn_x_4],[ttx,ratio]);
tx = double(soln.ttx);
throttle = double(soln.ratio);

% create time grid points
t = linspace(tmin,tmin+tMax,res);
u = nan(2,length(t));

% horizontal control

    u(1,t<(tx+tmin))  = xctrls(1) * throttle;
    u(1,t>=(tx+tmin)) = xctrls(2) * throttle;

% vertical control
    
    u(2,t<(ty+tmin))  = yctrls(1);
    u(2,t>=(ty+tmin)) = yctrls(2);
    
t = t + tMin;

disp(['x: ini: ' num2str(tMin) ...
      ' | mid: ' num2str(tx+tMin) ...
      ' | end: ' num2str(tMax+tMin) newline ...
      '   control: ' num2str(xctrls * throttle)])
disp(['y: ini: ' num2str(tMin) ...
      ' | mid: ' num2str(ty+tMin) ...
      ' | end: ' num2str(tMax+tMin) newline ...
      '   control: ' num2str(yctrls)])

end

