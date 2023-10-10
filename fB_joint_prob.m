function [states, prob_joint] = fB_joint_prob(data)

%% all possible states
n = size(data,2);  prob_integer = []; states = [];
for i = 1:2^n
    for j = 1:n
        comb(i, j)  = mod(floor((i-1)/2^(j-1)), 2);
    end
    
    ind = find(ismember(data, comb(end, :),'rows'));
   
    if ~isempty(ind)
    prob_integer = [prob_integer; numel(ind)];
    states = [states; unique(data(ind, :), 'rows')];
    end
end

prob_joint = prob_integer./size(data, 1);

