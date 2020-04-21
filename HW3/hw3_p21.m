close all
clear;clc

N = 6;

rw_map = -Inf * ones(7,7);
mv_map = repmat({' '},7,7);

s0 = [1 1]; % start point
rw_map(7,1) = 0; % end point

maps = repmat({' '},7,7,N+1);

for t = 0:N
    T = N-t;
    [rw_map,mv_map] = V(s0,T,rw_map,mv_map);
    maps(:,:,T+1) = mv_map;
end

row = s0(1);
col = s0(2);

reward = rw_map(row,col)
path = upper(showpath(maps,s0))


%% part (i)

actions = ['n','s','e','w','p'];

rewards_map = zeros(7,7,5);

for k = 1:5
    for i = 1:7
        for j = 1:7
            [~,rewards_map(i,j,k)] = do_action([i,j],actions(k));
        end
    end
    latex_table = latex(vpa(sym(rewards_map(:,:,k)),3))
end

% profile viewer

%% function definitions

function [next_state,reward] = do_action(state,dir)

    next_state = state;

    map = [ 0 1.3 0.2 0.2 2.1 1 0.6;
            5.5 4.5 3.9 7.6 5.2 5.7 1.3;
            7.2 6.4 5.3 0.2 8.4 6.4 7.1;
            0 2.2 0.8 1.2 9.8 8.8 7.7;
            9.4 8.9 6.2 1.9 8.2 4.4 6.5;
            9.0 7.2 0.8 0.4 3.0 0.4 3.5;
            0.4 0.6 0.4 4.9 0.3 1.0 2.4 ];
        
    switch dir
        case 'n'
            action = [-1 0];
        case 's'
            action = [1 0];
        case 'e'
            action = [0 1];
        case 'w'
            action = [0 -1];
        case 'p'
            action = [0 0];
    end

    row = state(1) + action(1);
    col = state(2) + action(2);
    
    if row > size(map,1) || row < 1
        reward = -Inf;
        return
    end
    if col > size(map,2) || col < 1
        reward = -Inf;
        return
    end
    
    reward = action(1)^2 + action(2)^2 + ...
             ( map(row,col) - map(state(1),state(2)) )^2;
    reward = -reward;
         
    next_state = [row,col];
        
end

function [rw_map_next,mv_map_next] = V(s0,T,rw_map,mv_map)

    % define actions
    actions = ['n','s','e','w','p'];
    
    % get row and column sizes to iterate through
    rows = size(rw_map,1);
    cols = size(rw_map,2);
    
    % copy the reward and movement map to output
    rw_map_next = rw_map;
    mv_map_next = mv_map;

    % iterate through all elements on the grid
    for m = 1:rows
        for n = 1:cols
            
            % iterate through all possible actions at the element
            rw_max = -Inf; % first, set the max reward to zero
            state = [m,n];
            for i = 1:length(actions)
                
                % if this is the first step, only s0 can contain
                % a reward more than -Inf
                if T == 0 && ~all(state==s0)
                    rw = -Inf;
                else
                    [next_state,rw] = do_action(state,actions(i));
                    row = next_state(1);
                    col = next_state(2);
                    rw = rw + rw_map(row,col);
                end
                if rw > rw_max
                    rw_max = rw;
                    rw_map_next(m,n) = rw;
                    mv_map_next{m,n} = actions(i);
                end
            end
            
        end
    end
    

end

function path = showpath(maps,s0)

    row = s0(1);
    col = s0(2);
    
    path = '';

    for i = 1:size(maps,3)-1
        
        next_state = do_action([row,col],maps{row,col,i});
        path = [path, maps{row,col,i}, ' '];
        row = next_state(1);
        col = next_state(2);
        
    end

end