function CMI = fB_MI_cond(states, prob_joint, pair, N_trains)

others = setdiff([1:N_trains], pair); 

%% conditional MI
MI_con = [];
for ii = 1 : size(states, 1) % explore all possible combinations, states_prob_joint already unique vector!!!!
    
    states_temp = states(ii, :);
    
    others_state_temp = states_temp(others);
    
    p_ind = find(ismember(states(:, others), others_state_temp, 'rows'));
    p_joint_others = sum(prob_joint(p_ind)); 
    
    p_joint = prob_joint(ii, :);
   
    p_joint_pair_others = [];
    for iii = 1: numel(pair)
        pair_others_state_temp = states_temp([pair(iii) others]);
        p_ind = find(ismember(states(:, [pair(iii) others]), pair_others_state_temp, 'rows'));
        p_joint_pair_others = [p_joint_pair_others; sum(prob_joint(p_ind))];
        
    end
    
    if (p_joint ~=0 && p_joint_others ~=0 && p_joint_pair_others(1) ~=0 && p_joint_pair_others(2) ~=0)
        MI_con = [MI_con; p_joint*log2((p_joint*p_joint_others)/(p_joint_pair_others(1)*p_joint_pair_others(2)))];
    end
    

end

CMI = sum(MI_con);