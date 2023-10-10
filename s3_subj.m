clear; close all; clc

% This script demonstates how to estimate the B-index rate and the corresponding
% B-structure applied on a real-data subject. The script as input parameter: 
% number of surrogates (num_surr), and  produces the B-index matrix (B) and 
% its corresponding matrix of the B-structure (Icol).

%% Input parameters
numsurr = 10;

%% load 5 time series
pmax = 12; q = 12;  alpha = 0.05; SelCrit = 'aic';

path = ['D:\01_Research\11_B_index\BIND toolbox\'];
alldata = load(strcat(path, 'T017.txt'));
sbp = alldata(:,1);
rr = alldata(:,2);
pvr = alldata(:,3);
co = alldata(:,4);
dbp = alldata(:,6);
data = [rr pvr co sbp dbp];

%% B-index estimation and B-structure
ret = lrp_BindexRate(data, pmax, q, numsurr, alpha, SelCrit);
B = ret.B; % B-index
Icol = ret.Icol; % B-structure
sB = tril(B')+triu(B);
sIcol = tril(Icol')+triu(Icol);


