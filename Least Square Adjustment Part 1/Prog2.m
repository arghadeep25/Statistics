%--------------------------------------------------------------------------
%   
%   ADJUSTMENT CALCULATION I
%   Author         : Arghadeep Mazumder
%   Task           : 2
%  
%--------------------------------------------------------------------------

 clc;
 clear all;
 close all;
 
 %--------------------------------------------------------------------------
 %   Observations and initial values for unknowns
 %--------------------------------------------------------------------------
 data=load('Homework3-2.txt');
 x=data(:,1);
 y=data(:,2);
  %Vector of Observations
 L = y;
 %Initial values for unknowns
 a_0=1;
 b_0=1;
X_0 =[a_0; b_0]
% %Number of observations
 no_n =length(L);
 % %Number of unknowns
 no_u =length(X_0);
%Redundancy
 r =no_n-no_u;
 %--------------------------------------------------------------------------
 %  stochastic model
 %--------------------------------------------------------------------------
 %VC Matrix of the observations
S_LL =eye(no_n);
%Theoretical standard deviations
 sigma_0 =1
%Cofactor matrix of the observations
 Q_LL =1/sigma_0^2*S_LL;
%Weight matrix
 P =inv(Q_LL);
% %--------------------------------------------------------------------------
% %  Adjustment
% %--------------------------------------------------------------------------
 %break-off conditions
  epsilon =10^-5;
  delta =10^-15;
  max_x_hat = Inf;
 max_Phi = Inf;
 % %Number of iterations
 iteration=0;
 %Iteration
 while max_x_hat>epsilon || max_Phi>delta
    plot(x,y,'.')
     hold on
     plot(x,a_0*exp(b_0*x),'r')
     xlabel('x')
ylabel('Y(x)')
     figure
    %Vector of reduced observations
    L_0 =a_0.*exp(b_0*x);
     l =L-L_0;
    %Designmatrix
   A=[exp(b_0*x) a_0*x.*exp(b_0*x)];
     %Normal matrix
    N =A'*P*A;
   %Vector of absolute values
   n =A'*P*l;
   %Inversion of normal matrix / Cofactor matrix of the unknowns
   Q_xx =inv(N);
  %Solution of normal equation
    x_hat = Q_xx*n;
  %Adjusted unknowns
   X_hat =X_0+x_hat
  %Update
  X_0 =X_hat;
   a_0=X_hat(1);
    b_0=X_hat(2);
%   Check 1
    max_x_hat =max(abs(x_hat));
  %Vector of residuals
     v =A*x_hat-l;
 %Vector of adjusted observations
    L_hat =L+v;
  %Check 2
    Phi=a_0*exp(b_0*x);
     max_Phi=max(abs(L_hat-Phi));
    %Update number of iterations
     iteration=iteration+1;
      end
 %Empirical reference standard deviation
s_0 = sqrt(v'*P*v/r);
%VC matrix of adjusted unknowns
Qxx_hat = inv(N);
S_XX_hat = s_0^2 * Qxx_hat;
%Standard deviation of the adjusted unknows
s_X = sqrt(diag(S_XX_hat))
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
plot(x,y,'.')
xlabel('x')
ylabel('Y(x)')
figure
bar(v)

