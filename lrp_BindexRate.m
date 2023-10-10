%% Computation of dynamic B-index rate (from MIR and cMIR) for all pairs of series in an input data matrix

function ret=lrp_BindexRate(data,pmax,q,numsurr,alpha,SelCrit)
% clear; close all; clc; 
% pmax=20; q=30; numsurr=100; alpha=0.05; SelCrit='aic';
% percorso='C:\Users\Luca\WORK\RESEARCH\DATA\Javorka-HUT_MA\used_in_mypapers\DATA\';
% cartella='\1_phase\';
% nomefile='1mj070.txt';
% data=load([percorso cartella nomefile]);
[N,M]=size(data);
X=data-mean(data);

[pottaic,pottbic]=mos_idMVAR(X',pmax,0);
switch SelCrit
    case 'aic'
        p=pottaic;
    case 'bic'
        p=pottbic;
end


out=lrp_idVAR(X,1:M,1:M,1:p);
Am=out.eA'; Su=out.es2u;
% [Am,Su]=idMVAR(X',p,0);

[oIxy,oIxy_z]=lrp_MIR_CMIR(Am,Su,q); % original MIR and CMIR


data_s=zeros(N,M);
for is=1:numsurr
    for m=1:M
        data_s(:,m) = surr_iaafft(data(:,m));
    end
    X_s=data_s-mean(data_s);
    outsurr=lrp_idVAR(X_s,1:M,1:M,1:p); %here I am keeping the original model order
    surrAm=outsurr.eA'; surrSu=outsurr.es2u;
    
    [Ixytmp,Ixy_ztmp]=lrp_MIR_CMIR(surrAm,surrSu,q);
    sIxy(:,:,is)=Ixytmp;
    sIxy_z(:,:,is)=Ixy_ztmp;
end
threshMI=prctile(sIxy,100*(1-alpha),3);
threshCMI=prctile(sIxy_z,100*(1-alpha),3);


Ixy=oIxy; Ixy(Ixy<threshMI)=0; %thresholded MI
Ixy_z=oIxy_z; Ixy_z(Ixy_z<threshCMI)=0; %thresholded CMI

Ixy_bin=Ixy>0; %binary MI
Ixy_z_bin=Ixy_z>0; %binary CMI

B=(Ixy-Ixy_z)./(max(Ixy,Ixy_z)); %B-index (computed on thresholded MI,CMI)

Ibin=(Ixy>0 & Ixy_z>0); %"binary" network structure
Icol=Ibin.*B; %"colored" network structure
Icol(Icol==0)=NaN;

ret.oIxy=oIxy; % original MI
ret.oIxy_z=oIxy_z; % original CMI
ret.Ixy=Ixy; % thresholded MI
ret.Ixy_z=Ixy_z; % thresholded MI
ret.Ixy_bin=Ixy_bin; % binary MI
ret.Ixy_z_bin=Ixy_z_bin; % binary CMI
ret.B=B; % B-index
ret.Ibin=Ibin; % "binary" network structure
ret.Icol=Icol; % "colored" network structure
ret.p=p; %model order