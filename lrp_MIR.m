%% MIR from linear regression
% performs linear regression of the present state of N target processes from the past states of M driver processes
% i: indexes of the first process, X
% j: indexes of the second process, Y
% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)

function ret = lrp_MIR(Am,Su,q,i,j)

% SigmaX_XY=Su(i,i); SigmaY_XY=Su(j,j); %full model

[SigmaY,SigmaY_Y] = lrp_LinReg(Am,Su,q,j,j); % restricted model, Ypresent given Ypast
[SigmaX,SigmaX_X] = lrp_LinReg(Am,Su,q,i,i); % restricted model, Xpresent given Xpast
[SigmaXY,SigmaXY_XY] = lrp_LinReg(Am,Su,q,[i j],[i j]); % restricted model, XYpresent given XYpast

Hx=0.5*log(((2*pi*exp(1))^size(SigmaX,1))*det(SigmaX_X)); %entropy rate of X
Hy=0.5*log(((2*pi*exp(1))^size(SigmaY,1))*det(SigmaY_Y)); %entropy rate of Y
Hxy=0.5*log(((2*pi*exp(1))^size(SigmaXY,1))*det(SigmaXY_XY)); %entropy rate of [X Y]
Ixy=Hx+Hy-Hxy; % MIR

ret.Hx=Hx;
ret.Hy=Hy;
ret.Hxy=Hxy;
ret.Ixy=Ixy;

