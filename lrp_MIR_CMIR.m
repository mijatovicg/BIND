%% Computation of dynamic MIR and CMIR (linear Gaussian) for all pairs of series in a multivariate process
%% given VAR parameters, computes MIR between each pair of processes, and conditional MIR given all other processes
% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)


function [Ixy,Ixy_z]=lrp_MIR_CMIR(Am,Su,q)
% clear; close all; clc;
% Q=4;
% par.poles{1}=([0.8 0.3]);
% par.poles{2}=([0.8 0; 0.82 0.1]);
% par.poles{3}=[0.9 0.1];
% par.poles{4}=[];
% par.Su=[4 2 3 2];
% par.coup=[1 2 1 0.8; 1 3 1 0.4; 2 3 2 0.6; 2 4 3 -0.7];
% q=30;
% [Am,Su,Ak]=theoreticalVAR(Q,par); % % VAR parameters
% Am=Am';

Q=size(Su,1); %n. of processes

Ixy=zeros(Q,Q); Ixy_z=Ixy; %initialize MIR and cMIR
for ix=1:Q
    for iy=ix+1:Q
        iz=setdiff(1:Q,[ix iy]);
        retx_y=lrp_MIR(Am,Su,q,ix,iy);
        Ixy(ix,iy)=retx_y.Ixy; %MIR btw processes ix and iy
        
        retx_z=lrp_MIR(Am,Su,q,ix,iz); %retx_z.Ixy; %MIR btw processes ix and iz
        retx_yz=lrp_MIR(Am,Su,q,ix,[iy iz]); %retx_yz.Ixy; %MIR btw processes ix and [iy,iz]
        Ixy_z(ix,iy)=retx_yz.Ixy-retx_z.Ixy; %cMIR btw processes ix and iy given processes iz  
    end
end