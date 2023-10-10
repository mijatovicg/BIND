function MI = fB_MI(states, prob_joint, pair)

subset = states(:, pair);
subset_un = unique(subset, 'rows');

MI_array = [];
for iii = 1: size(subset_un, 1)
    subset_un_tmp = subset_un(iii, :);
    p_ind = find(ismember(subset, subset_un_tmp, 'rows'));
    
    p_joint_subset = sum(prob_joint(p_ind)); % joint probability of the subset
    
    p_marg_all = [];
    for iiii = 1: numel(subset_un_tmp)
        p_ind = find(subset(:, iiii) == subset_un_tmp(iiii));
        p_marg_all = [p_marg_all sum(prob_joint(p_ind))]; % marginal probability of the target
    end
    
    if (p_joint_subset ~=0 && p_marg_all(1) ~=0 && p_marg_all(2) ~=0)
        MI_array = [MI_array; p_joint_subset*(log2((p_joint_subset /(p_marg_all(1)*p_marg_all(2)))))];
    end
    
end

MI = sum(MI_array);


