%--------------------------------------------------------------------------
%   
%   ADJUSTMENT CALCULATION I
%   Author         : Arghadeep Mazumder
%   Task           : 1
%  
%--------------------------------------------------------------------------
clc;
clear all;
close all;
%--------------------------------------------------------------------------
%             Observations and Initial Values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
data=load('Homework3-1.txt');
t=data(:,1); % Time axis
y=data(:,2); %Converted to radian

% %Vector of observations
L = y;
% number of observation
no_n=length(L);
% %Number of unknowns
 no_u = 5;
% 
% %Redundancy
 r = no_n-no_u;
 %--------------------------------------------------------------------------
%                              Stochastic Model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL =eye(no_n);
%Theoretical standard deviations
 sigma_0 =1;
%Cofactor matrix of the observations
 Q_LL =1/sigma_0^2*S_LL;
%Weight matrix
 P =inv(Q_LL);
 %--------------------------------------------------------------------------
 %                                 Adjustment
 %--------------------------------------------------------------------------
 %Designmatrix
A = [t.^4 t.^3 t.^2 t ones(no_n,1)]
% %Normal matrix
N =A'*P*A;
% %Vector of absolute values
 n =A'*P*L;
% %Inversion of normal matrix / Cofactor matrix of the unknowns
 Q_xx_hat = inv(N);
% %Adjusted unknowns
 X_hat = Q_xx_hat*n;
Unknowns=X_hat 
% %Vector of residuals
 v = A*X_hat-L;
% %Vector of adjusted observations
 L_hat = L + v;
%Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);
%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx_hat;
%Standard deviation of the adjusted unknows
Standard_Deviation = sqrt(diag(S_XX_hat))
%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx_hat*A';
%VC matrix of adjusted observations
S_LL_hat =s_0^2*Q_LL_hat;
%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;
%VC matrix of residuals
S_vv = s_0^2*Q_vv;
%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
 plot(t,L_hat,'.',t,y,'.')
  xlabel('time')
 ylabel('Inclination')
 figure
 plot(v)
 xlabel('Observstion')
 ylabel('Residual')
 

