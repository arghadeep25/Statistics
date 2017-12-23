clc;
clear all;
close all;
format long g

% Measured 2D coordinates
set1=load('Straightline_set1.txt');
set2=load('Straightline_set2.txt');
y=[set1(:,2);set2(:,2)];
x=[set1(:,1);set2(:,1)];
L = [y];

%Number of condition equations
no_n=length(L);

%Number of unknowns
no_u=2;

%Redundancy
r=no_n-no_u;

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%Weight matrix
S_LL= eye(length(L));
sigma_0 = 1;
Q_LL= 1/sigma_0^2 *S_LL;
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Designmatrix A
A = [x,ones(no_n,1)];

%Normal matrix
N= A'*P*A;

%Vector of absolute values
n= A'*P*L;
    
% Inversion of normal matrix / Cofactor matrix of the unknowns
Q_xx=inv(N);
    
%Adjusted unknowns
X_hat=Q_xx*n;
        
% Vector of residuals
v=A*X_hat-L;

% Vector of adjusted observations
L_hat=L+v;

% Final Check
if (L_hat-(X_hat(1)*x+X_hat(2)))<10^-15
    disp('Everything fine')
else
    disp('There is a problem')
end

% Empirical reference standard deviation
s_0=sqrt(v'*P*v/r);

% VC matrix of adjusted unknowns
S_XX_hat=(s_0^2)*Q_xx;

% Standard deviation of the adjusted unknowns
s_X=sqrt(diag(S_XX_hat));

% Cofactor matrix of adjusted observations
Q_LL_hat=A*Q_xx*A';

% VC matrix of adjusted observations
S_LL_hat=(s_0^2)*Q_LL_hat;

% Standard deviation of the adjusted observations
s_L_hat=sqrt(diag(S_LL_hat));

% Cofactor matrix of the residuals
Q_vv=Q_LL-Q_LL_hat;

% VC matrix of residuals
S_vv=(s_0^2)*Q_vv;

% Standard deviation of the residuals
s_v=sqrt(diag(S_vv));

figure
bar(v*1000)
xlabel('Number of Residual')
ylabel('v (mm)')
title('Residuals')

xseries=min(x):0.1:max(x);
figure
plot(set1(:,1),set1(:,2),'ro'), hold on
plot(set2(:,1),set2(:,2),'b*'), hold on
plot(xseries,xseries*X_hat(1)+X_hat(2),'k'), hold on
xlabel('X Direction [m]')
ylabel('Y DIrection [m]')
title('Data Points & Adjusted Line')
legend('SET 1 Points','SET 2 Points','ADJUSTED LINE','location','southeast')
