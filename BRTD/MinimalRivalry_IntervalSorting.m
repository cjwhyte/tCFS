% Interval sorting script

function [X_dom_FR, X_sup_FR, X_dom_probe_FR, X_sup_probe_FR] = MinimalRivalry_IntervalSorting(X_store,p)
   
    % n sorting
    num_intervals = 6;
    summation_length = 10/p.DT;

    % grab switch times
    X_switch_times = find(abs(diff(X_store(3,:)))>0);
    
    % storage structures
    X_struct = {}; 
    
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
            X_struct{counter} = X_store(:,start_idx:end_idx);
        end 
    end 
   
    dom_counter = 0; sup_counter = 0; 
    dom_probe_counter = 0; sup_probe_counter = 0; 
    % sort activity into invervals of 6
    for ii = 1:size(X_struct,2)
        trl_length = size(X_struct{ii},2);
        interval_size =  floor(trl_length/num_intervals);
        intervals = 1:interval_size:interval_size*num_intervals+1;
        % grab activity when R is dominant and no probe
        if sum(X_struct{ii}(3,:))>0 && isempty(find(X_struct{ii}(4,:)>0))
            dom_counter = dom_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(X_struct{ii},2)
                    end_idx = size(X_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                rand_time_idx = randi([intervals(1),end_idx-summation_length]);
                X_dom_FR(dom_counter,jj) = sum(X_struct{ii}(2,rand_time_idx:rand_time_idx+summation_length)).*p.DT;
            end 
        % grab activity when R is suppressed and no probe
        elseif sum(X_struct{ii}(3,:))<0 && isempty(find(X_struct{ii}(4,:)>0))
            sup_counter = sup_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(X_struct{ii},2)
                    end_idx = size(X_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                rand_time_idx = randi([intervals(1),end_idx-summation_length]);
                X_sup_FR(sup_counter,jj) = sum(X_struct{ii}(2,rand_time_idx:rand_time_idx+summation_length)).*p.DT;
            end 
        % grab activity when R is dominant and there is a probe
        elseif sum(X_struct{ii}(3,:))>0 & ~isempty(find(X_struct{ii}(4,:)>0)) & find(X_struct{ii}(4,:)>0) + summation_length < size(X_struct{ii}(4,:),2)
            dom_probe_counter = dom_probe_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(X_struct{ii},2)
                    end_idx = size(X_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                X_dom_probe_norm{dom_probe_counter,jj} = X_struct{ii}(2,intervals(jj):end_idx);
            end 
            probe_idx = find(X_struct{ii}(4,:)>0);
            X_dom_probe_FR(dom_probe_counter,1) = sum(X_struct{ii}(2,probe_idx(1):probe_idx(1)+summation_length)).*p.DT;
            [~,X_dom_probe_FR(dom_probe_counter,2)] = min(abs(intervals(2:end)-probe_idx(1)),[],"All");
        % grab activity when R is suppressed and there is a probe
        elseif sum(X_struct{ii}(3,:))<0 & ~isempty(find(X_struct{ii}(4,:)>0)) & find(X_struct{ii}(4,:)>0) + summation_length < size(X_struct{ii}(4,:),2)
            sup_probe_counter = sup_probe_counter + 1;
            for jj = 1:num_intervals
                if intervals(jj+1) > size(X_struct{ii},2)
                    end_idx = size(X_struct{ii},2);
                else 
                    end_idx = intervals(jj+1);
                end 
                X_sup_probe_norm{sup_probe_counter,jj} = X_struct{ii}(2,intervals(jj):end_idx);
            end 
            probe_idx = find(X_struct{ii}(4,:)>0);
            X_sup_probe_FR(sup_probe_counter,1) = sum(X_struct{ii}(2,probe_idx(1):probe_idx(1)+summation_length)).*p.DT;
            [~,X_sup_probe_FR(sup_probe_counter,2)] = min(abs(intervals(2:end)-probe_idx(1)),[],"All");
        end 
    end 

end 