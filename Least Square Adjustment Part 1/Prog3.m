%--------------------------------------------------------------------------
%   
%   ADJUSTMENT CALCULATION I
%   Author         : Arghadeep Mazumder
%   Task           : 3
%  
%--------------------------------------------------------------------------
clc;
clear all;
close all;
%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
L = [ 1.160; 15.15];
x=8.93;
%Initial values for unknowns
V_0=1.6;
%Number of observations
no_n = length(L);
%Number of unknowns
no_u = 1;
%Redundancy
r = no_n-no_u;
%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL =[0.005^2 0.0 ;0.0 0.05^2 ];
%Theoretical standard deviation
sigma_0 =1;
%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;
%Weight matrix
P = inv(Q_LL);
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-15;
max_x_hat =Inf;
max_Phi= Inf;
%Number of iterations
iteration=0;
%Iteration
while max_x_hat>epsilon || max_Phi > delta

    %Vector of reduced observations
    L_0(1,1)=(V_0)^(1/3);
    L_0(2,1)= x*V_0;
    l = L-L_0;
   %Designmatrix
    A=[1/3*(V_0)^(-2/3); x];
   %Normal matrix
    N =A.'*P*A;
    %Vector of absolute values
    n =A'*P*l;
    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx =inv(N);
    %Solution of ;normal equation
    x_hat =Q_xx*n;
    %Adjusted unknowns
    X_hat =V_0+x_hat;
   %Update
    V_0 =X_hat;
   %Check 1
    max_x_hat =max(abs(x_hat));
    %Vector of residuals
    v =A*x_hat-l;
    %Vector of adjusted observations
    L_hat =L+v;
    %Check 2
    Phi(1,1)=(V_0)^(1/3);
    Phi(2,1)= x*V_0;
    max_Phi =max(abs(L_hat-Phi));
    %Update number of iterations
    iteration=iteration+1;
end
Volume= V_0
 %Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);
%VC matrix of adjusted unknowns
Qxx_hat = inv(N);
S_XX_hat = s_0^2 * Qxx_hat;
%Standard deviation of the adjusted unknows
Standard_Deviation = sqrt(diag(S_XX_hat))
%Cofactor matrix of adjusted observations
Q_LL_hat = A *Qxx_hat *A';
%VC matrix of adjusted observations
S_LL_hat =s_0^2 * Q_LL_hat;
%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat,0));
%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;
%VC matrix of residuals
S_vv = sigma_0^2 * Q_vv;
%Standard deviation of the residuals
s_v = diag(S_vv,0);




    
