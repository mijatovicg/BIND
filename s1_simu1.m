clear; close all; clc;

% This script generates a single realization of the simulated dataset constited of binary 
% variables based on the simulation description given in the section B of the article
% "Assessing High-Order Links in Cardiovascular and Respiratory Networks via Static and Dynamic Information Measures" 
% by G. Mijatovic, L. Sparacino, Y. Antonacci, M. Javorka, D. Marinazzo, S. Stramaglia, and L. Faes.

% It takes as input parameters: 
% number of observations for the simulated dataset (N), 
% number of surrogates (num_surr), 
% probabilities (gamma1, gamma2, and gamma3) which modulate the strength of coupling. 

% The script produces:
% the B-index matrix (B) and 
% its corresponding matrix of the B-structure (Icol).

%% input parameters
N = 100; % number of observations
num_surr = 100;% number of surrogates
gamma1 = 0.9; 
gamma2 = 0.9; 
gamma3 = 0.8; % probabilities
%% generate a single realization of a dataset
X1 = randi([0,1], N, 1); X3 = randi([0,1], N, 1); X4 = randi([0,1], N, 1); X5 = randi([0,1], N, 1); X9 = randi([0,1], N, 1); % i.i.d. variables

X2 = [];
for ii = 1 : N
    if (rand(1) <= gamma1 && X3(ii)+ X4(ii) + X5(ii) == 0)
        X2 = [X2; 0];
    else
        X2 = [X2; 1]; % 3-inputs OR gate
    end
end

n_ind = round(numel(X5) - numel(X5)*gamma2);
r = randi([1 numel(X5)],1, n_ind);
X6 = X5;
X6(r) = not(X6(r)); % noisy copy

n_ind = round(numel(X5) - numel(X5)*gamma2);
r = randi([1 numel(X5)],1, n_ind);
X7 = X5;
X7(r) = not(X7(r)); % noisy copy

X8 = [];
for ii = 1 : N
    if (rand(1) <= gamma1)
        X8 = [X8; or(X6(ii), X7(ii))];
    else
        X8 = [X8; not(or(X6(ii), X7(ii)))]; %2-inputs OR gate
    end
end

n_ind = round(numel(X9) - numel(X9)*gamma3);
r = randi([1 numel(X9)],1, n_ind);
X10 = X9;
X10(r) = not(X10(r)); % noisy copy

data = [X1 X2 X3 X4 X5 X6 X7 X8 X9 X10]; % collect all variables in data
N_trains = size(data,2);

%% B-index estimation and B-structure
B = nan(N_trains); sIcol = B;
for i = 1 : size(data, 2)
    
    for j = i+1 : size(data, 2)
        
        pair = [i, j];
        [states, prob_joint] = fB_joint_prob(data);
        MI = fB_MI(states, prob_joint, pair);
        CMI = fB_CMI(states, prob_joint, pair, N_trains);
        
        %% surrogates based on random shuffling
        cnt = 1; MI_surr = []; CMI_surr = [];
        while cnt <= num_surr
            
            tmp_i = data(:, i);
            tmp_j = data(:, j);
            
            data_perm = data;
            data_perm(:, i) = tmp_i(randperm(numel(tmp_i)));
            data_perm(:, j) = tmp_j(randperm(numel(tmp_j)));
            
            [states_perm, prob_joint_perm] = fB_joint_prob(data_perm);
            MI_surr = [MI_surr; fB_MI(states_perm, prob_joint_perm, pair)];
            CMI_surr = [CMI_surr; fB_CMI(states_perm, prob_joint_perm, pair, N_trains)];
            
            cnt = cnt+1;
        end
       
        thr_MI = prctile(MI_surr, 95);
        thr_CMI = prctile(CMI_surr, 95);
      
        %%
        if MI < thr_MI 
            MI = 0;
        end
        if CMI < thr_CMI
            CMI = 0;
        end
        
        B_val = (MI - CMI)/max(MI, CMI);% B index
        B(i, j) = B_val; B(j, i) = B_val;
        
        if B_val == 1
            sIcol(i, j) = NaN;  sIcol(j, i) = NaN;
        elseif B_val == -1
            sIcol(i, j) = NaN; sIcol(j, i) = NaN;
        else
            sIcol(i, j) = B_val;  sIcol(j, i) = B_val; % B_structure based on the B-index
        end
        
    end % jth train
end % ith train



