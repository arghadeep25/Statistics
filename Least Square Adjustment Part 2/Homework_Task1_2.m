clc
close all
clear all
%Vector of Observations
L = [ 10.509; 5.36; -8.523; -7.348; -3.167; 15.881];
%Fixed
H_B=455.873;
%Observation Vector
benchmarks = [H_B; -H_B; 0; 0; -H_B;0]; 
L_dash = L - benchmarks;
%Number of observations
no_n = 6;
%Number of unknowns
no_u =3;
%Redundancy
r = no_n - no_u; 
%Cofactor matrix of the observations
S_LL =diag([0.006 0.004 0.005 0.003 0.004 0.012]).^2;
%Theoretical standard deviation
sigma_0 =1;
%Cofactor matrix of the observations
Q_LL =1/sigma_0^2*S_LL;
%Weight matrix
P =inv(Q_LL);
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Designmatrix
A = [-1  0 0 ; 
    0 1 0; 
    0 -1 1;   
    1 0 -1;    
    0 0 1;    
    -1 1 0];
%Normal matrix
N = A'*P*A;
%Vector of absolute values
n = A'*P*L_dash;
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_xx = inv(N);
%Adjusted unknowns
X_hat =Q_xx*n;
%Vector of residuals
v=A*X_hat-L_dash;
%Vector of adjusted observations
L_hat=L+v;
%Check 
delta =10^-10;
if max(abs(L_hat-(A*X_hat+benchmarks)))<delta
    disp('everything is fine')
else
    disp('We have a Problem')
end
Adjusted_Unknown=X_hat
Residuals=v
Adjusted_Observation=L_hat
%Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);
%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;
%Standard deviation of the adjusted unknows
Std_Dev_Adj_Unkwn = sqrt(diag(S_XX_hat))
%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';
%VC matrix of adjusted observations
S_LL_hat =s_0^2*Q_LL_hat;
%Standard deviation of the adjusted observations
Std_Dev_Adj_Obs = sqrt(diag(S_LL_hat))
%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;
%VC matrix of residuals
S_vv = s_0^2*Q_vv;
%Standard deviation of the residuals
Std_Dev_Residuals = sqrt(diag(S_vv))