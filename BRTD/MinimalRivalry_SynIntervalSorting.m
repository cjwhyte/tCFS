% Interval sorting script

function [SD_dom, SD_sup] = MinimalRivalry_SynIntervalSorting(X_store,SD_store,p)
   
    % n sorting
    num_intervals = 6;

    % grab switch times
    X_switch_times = find(abs(diff(X_store(3,:)))>0);
    
    % storage structures
    X_struct = {}; SD_struct = {}; 
    
    % sort dominance intervals
    counter = 0;
    for ii = 1:length(X_switch_times)-1
        if ii == 1
           start_idx = 1;
           end_idx = X_switch_times(ii);
        else 
           start_idx = X_switch_times(ii)+1;
           end_idx = X_switch_times(ii+1);
        end 
        if (end_idx-start_idx) > 600/p.DT
            counter = counter + 1;
            SD_struct{counter} = SD_store(start_idx:end_idx);
            X_struct{counter} = X_store(:,start_idx:end_idx);
        end 
    end 
   
    dom_counter = 0; sup_counter = 0; 
    % sort activity into invervals of 6
    for ii = 1:size(X_struct,2)
        trl_length = size(X_struct{ii},2);
        interval_size =  floor(trl_length/num_intervals);
        intervals = 1:interval_size:interval_size*num_intervals+1;
        % grab suppression depth when R is dominant and no probe
        if sum(X_struct{ii}(3,:))>0
            dom_counter = dom_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(SD_struct{ii},2)
                    end_idx = size(SD_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                SD_dom(dom_counter,jj) = mean(SD_struct{ii}(intervals(jj):end_idx));
            end 
        % grab suppression depth when R is suppressed 
        elseif sum(X_struct{ii}(3,:))<0
            sup_counter = sup_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(SD_struct{ii},2)
                    end_idx = size(SD_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                SD_sup(sup_counter,jj) = mean(SD_struct{ii}(intervals(jj):end_idx));
            end 
        end 
    end 

end 