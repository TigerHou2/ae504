close all
clear;clc

N = 3;

rw_map = -Inf * ones(5,1); % rewards at each state
mv_map = zeros(5,1); % action at each state

s0 = 1; % start point
rw_map(2,1) = 0; % end point

rw_maps = zeros(5,N+1);
mv_maps = zeros(5,N+1); % tracks actions at each state for each step

for t = 0:N
    T = N-t;
    [rw_map,mv_map] = V(s0,T,rw_map,mv_map);
    rw_maps(:,T+1) = rw_map;
    mv_maps(:,T+1) = mv_map;
end

[path,states,reward] = showpath(mv_maps,s0)


% profile viewer

%% function definitions

function [next_state,reward] = do_action(state,action)

    T = [ 3 1 5 2; ...
          2 3 4 4; ...
          2 5 4 6; ...
          2 7 4 8; ...
          2 9 4 10 ];
      
    next_state = T(state,2*action-1);
    reward = T(state,2*action);
        
end

function [rw_map_next,mv_map_next] = V(s0,T,rw_map,mv_map)

    % define actions
    actions = [1,2];
    
    % get all states to iterate through
    states = length(rw_map);
    
    % copy the reward and movement map to output
    rw_map_next = rw_map;
    mv_map_next = mv_map;

    % iterate through all elements
    for m = 1:states
            
            % iterate through all possible actions at the element
            state = m;
            for i = 1:length(actions)
                
                % if this is the first step, only s0 can contain
                % a reward more than -Inf
                if T == 1 && state~=s0
                    rw = -Inf;
                else
                    [next_state,rw] = do_action(m,actions(i));
                    rw = rw + rw_map(next_state);
                end
                if rw > rw_map_next(m)
                    rw_map_next(m) = rw;
                    mv_map_next(m) = actions(i);
                end
            end
            
    end
    

end

function [path,states,reward] = showpath(maps,s0)

    state = s0;
    
    path = [];
    states = [];
    reward = 0;

    for i = 2:size(maps,2)
        
        [next_state,rw] = do_action(state,maps(state,i));
        states = [states; state];
        path = [path; maps(state,i)];
        state = next_state;
        reward = reward + rw;
        
    end
    
    states = [states; do_action(states(end),maps(states(end),i))];

end