clear; close all; clc;

% This script generates a single realization of the simulated dataset constited of Gaussian 
% variables based on the simulation description given in the section B of the article
% "Assessing High-Order Links in Cardiovascular and Respiratory Networks via Static and Dynamic Information Measures" 
% by G. Mijatovic, L. Sparacino, Y. Antonacci, M. Javorka, D. Marinazzo, S. Stramaglia, and L. Faes.

% It takes as input parameters: 
% number of observations for the simulated dataset (N), 
% number of surrogates (num_surr), 
% time lagged coefficients (a_ij) 
% and the maximum imposed lag (p). 

% The script produces:
% the B-index matrix (B) and 
% its corresponding matrix of the B-structure (Icol).

%% Input parameters
N = 250; % number of observations
numsurr = 100; % number of surrogates
pmax = 2; % maximum lag
a_ij = 0; % time lagged coefficients
%%
Q = 6; % n. of processes
alpha = 0.05;  SelCrit = 'aic';
q = 20;
scen_array = [1, 2]; % leaves behave as: 1 – receivers, 2 - mediators
for sca = 1 : numel(scen_array)
    
    scen = scen_array(sca);
    
    if scen == 2 % propagating
        a21 = a_ij; a31 = a_ij; a41 = a_ij; a51 = a_ij;
        a62 = 1-a_ij; a63 = 1-a_ij; a64 = 1-a_ij; a65 = 1-a_ij;
        a26 = 0; a36 = 0; a46 = 0; a56 = 0;
    elseif scen == 1 % competing
        a21 = a_ij; a31 = a_ij; a41 = a_ij; a51 = a_ij;
        a62 = 0; a63 = 0; a64 = 0; a65 = 0;
        a26 = 1-a_ij; a36 = 1-a_ij; a46 = 1-a_ij; a56 = 1-a_ij;
    else
    end
   
    %% simulation design
    Su = eye(Q);
    Ak = zeros(Q, Q, pmax); % blocks of coefficients
    % effects originating from 1 (at lag 1)
    Ak(2,1,1) = a21;
    Ak(3,1,1) = a31;
    Ak(4,1,1) = a41;
    Ak(5,1,1) = a51;
    % effects directed to 6 (at lag 1)
    Ak(6,2,1) = a62;
    Ak(6,3,1) = a63;
    Ak(6,4,1) = a64;
    Ak(6,5,1) = a65;
    % effects originating from 6 (at lag 2)
    Ak(2,6,2) = a26;
    Ak(3,6,2) = a36;
    Ak(4,6,2) = a46;
    Ak(5,6,2) = a56;
    
    Am=[];
    for kk = 1 : pmax
        Am = [Am Ak(:,:,kk)];
    end
    
    %% B-index estimation and B-structure
    U = mvnrnd(zeros(Q,1), Su, N)'; % U: Q*N matrix of innovations
    Y = MVARfilter(Am, U); % simulated time series
    data = Y';
    
    ret = lrp_BindexRate(data, pmax, q, numsurr, alpha, SelCrit);
    B = ret.B; % B-index
    Icol = ret.Icol; % B-structure
    sB = tril(B')+triu(B);
    sIcol = tril(Icol')+triu(Icol); 
    
end % sca, propagation Vs. competition








