close all
clear;clc

s0 = 1;
sN = 2;
N = 10;

reward = V(s0,N,sN)

%% function definitions

function [states,actions] = load_sa

    states = [1 2 3 4 5];
    actions = [1 2];

end

function [next_state,reward] = do_action(state,action)

    switch state
        
        case 1
            if action == 1
                next_state = 3;
                reward = 1;
            elseif action == 2
                next_state = 5;
                reward = 2;
            end
            
        case 2
            if action == 1
                next_state = 2;
                reward = 3;
            elseif action == 2
                next_state = 4;
                reward = 4;
            end
            
        case 3
            if action == 1
                next_state = 2;
                reward = 5;
            elseif action == 2
                next_state = 4;
                reward = 6;
            end
            
        case 4
            if action == 1
                next_state = 2;
                reward = 7;
            elseif action == 2
                next_state = 4;
                reward = 8;
            end
            
        case 5
            if action == 1
                next_state = 2;
                reward = 9;
            elseif action == 2
                next_state = 4;
                reward = 10;
            end
            
    end

end

function reward = V(state,T,sN)
    
    [~,actions] = load_sa;
    
    reward = 0;

    if T == 0
        if state ~= sN
            reward = -Inf;
            disp([newline 'Bad path!' ...
                  ' (s=' num2str(state) ', T=' num2str(T) ')'])
        else
            disp([newline 'Good path!' ...
                  ' (s=' num2str(state) ', T=' num2str(T) ')'])
        end
        return
    end

    r = 0;
    rw = 0;
    for i = 1:length(actions)
        [~,ww] = do_action(state,actions(i));
        rr = ww + V(do_action(state,actions(i)),T-1,sN);
        if rr > r
            r  = rr;
            rw = ww;
        end
    end
    
    reward = r;
    
    disp(['T = ' num2str(T) ...
          ', state = ' num2str(state) ...
          ', a = ' num2str(actions(i)) ...
          ', reward = ' num2str(rw)])

end