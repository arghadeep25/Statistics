%--------------------------------------------------------------------------
%   
%   ADJUSTMENT CALCULATION I
%   Template for non-linear functional models
% 
%   Author         : Sven Weisbrich
%   Version        : April 13, 2012
%   Last changes   : January 29, 2014
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long;
%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations
dir = load('Direction_Observations.txt');
dist = load('Distance_Observations.txt');
dir = dir*pi/200;

L = [dist(:,3); dir(:,3)];

%Fixed Coordinates
y6=5317651.428;
x6=4968940.373;
y9=5324162.853;	
x9=4970922.160;

%Initial values for unknowns
y1=5314698.13;	
x1=4965804.18;
y15=5320448.85;	
x15=4962997.53;
w1=0;
w6=0;
w9=0;
w15=0;
X_0 =[x1 y1 x15 y15 w1 w6 w9 w15]'; 

    %Number of observations
    no_n = length(L);

    %Number of unknowns
    no_u = length(X_0);

    %Redundancy
    r = no_n-no_u;

%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sigma_s=0.1;
sigma_r=0.001*pi/200;
S_LL = diag([sigma_s^2*ones(5,1); sigma_r^2*ones(9,1)]);

%Theoretical standard deviation
sigma_0 =0.001;

    %Cofactor matrix of the observations
    Q_LL = 1/sigma_0^2*S_LL;

    %Weight matrix
    P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
    %break-off conditions
    epsilon=10^-5;
    delta=10^-9;
    max_x_hat=Inf;
    check2=Inf;

    %Number of iterations
    iteration=0;

%Iteration
while max_x_hat>epsilon || check2>delta

%Vector of reduced observations
L_0(1) = calcdist (x6, y6, x1, y1);
L_0(2) = calcdist (x9, y9, x1, y1);
L_0(3) = calcdist (x9, y9, x6, y6);
L_0(4) = calcdist (x15, y15, x1, y1);
L_0(5) = calcdist (x15, y15, x9, y9);

 %direction ik =direction(yk, xk, yi, xi, wi)
L_0(6) = direction (y6, x6, y1, x1, w1);
L_0(7) = direction (y15, x15, y1, x1, w1);    
L_0(8) = direction (y1, x1, y6, x6, w6);     
L_0(9) = direction (y9, x9, y6, x6, w6);     
L_0(10) = direction (y15, x15, y9, x9, w9);     
L_0(11) = direction (y1, x1, y9, x9, w9);     
L_0(12) = direction (y6, x6, y9, x9, w9);     
L_0(13) = direction (y1, x1, y15, x15, w15);     
L_0(14) = direction (y9, x9, y15, x15, w15); 

    l=L-L_0';

    %Designmatrix
    
A(1,1) = ds_dx_target(x6,y6,x1,y1);
A(1,2) = ds_dy_target(x6,y6,x1,y1);
A(2,1) = ds_dx_target(x9,y9,x1,y1);
A(2,2) = ds_dy_target(x9,y9,x1,y1);
A(4,1) = ds_dx_target(x15,y15,x1,y1);
A(4,2) = ds_dy_target(x15,y15,x1,y1);
A(4,3) = ds_dx_origin(x15,y15,x1,y1);
A(4,4) = ds_dy_origin(x15,y15,x1,y1);
A(5,3) = ds_dx_origin(x15,y15,x9,y9);
A(5,4) = ds_dy_origin(x15,y15,x9,y9);

%Directions
A(6,1) = dr_dx_origin(x1,y1,x6,y6);
A(6,2) = dr_dy_origin(x1,y1,x6,y6);
A(6,5) = -1;

A(7,1) = dr_dx_origin(x1,y1,x15,y15);
A(7,2) = dr_dy_origin(x1,y1,x15,y15);
A(7,3) = dr_dx_target(x1,y1,x15,y15);
A(7,4) = dr_dy_target(x1,y1,x15,y15);
A(7,5) = -1;

A(8,1) = dr_dx_target(x6,y6,x1,y1);
A(8,2) = dr_dy_target(x6,y6,x1,y1);
A(8,6) = -1;

A(9,6)=-1;

A(10,3) = dr_dx_target(x9,y9,x15,y15);
A(10,4) = dr_dy_target(x9,y9,x15,y15);
A(10,7) = -1;

A(11,1) = dr_dx_target(x9,y9,x1,y1);
A(11,2) = dr_dy_target(x9,y9,x1,y1);
A(11,7) = -1;

A(12,7) = -1;

A(13,1) = dr_dx_target(x15,y15,x1,y1);
A(13,2) = dr_dy_target(x15,y15,x1,y1);
A(13,3) = dr_dx_origin(x15,y15,x1,y1);
A(13,4) = dr_dy_origin(x15,y15,x1,y1);
A(13,8) = -1;

A(14,3) = dr_dx_origin(x15,y15,x9,y9);
A(14,4) = dr_dy_origin(x15,y15,x9,y9);
A(14,8) = -1;
    
    %Normal matrix
    N=A'*P*A;

    %Vector of absolute values
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx =inv(N);
    
    %Solution of normal equation
    x_hat = Q_xx*n;
    
    %Adjusted unknowns
    X_hat = X_0+x_hat;
    X_0=X_hat;
    
    %Update
x1=X_0(1);
y1=X_0(2);
x15=X_0(3);
y15=X_0(4);
w1=X_0(5);
w6=X_0(6);
w9=X_0(7);
w15=X_0(8);
    %Check 1
    max_x_hat=max(abs(x_hat));
    %Update number of iterations
    iteration=iteration+1;
    %Vector of residuals
    v = A*x_hat-l;
    %Vector of adjusted observations
    L_hat = L+v;

    %Check 2
Psi(1) = calcdist (x6, y6, x1, y1);
Psi(2) = calcdist (x9, y9, x1, y1);
Psi(3) = calcdist (x9, y9, x6, y6);
Psi(4) = calcdist (x15, y15, x1, y1);
Psi(5) = calcdist (x15, y15, x9, y9);

Psi(6) = direction (y6, x6, y1, x1, w1);
Psi(7) = direction (y15, x15, y1, x1, w1);    
Psi(8) = direction (y1, x1, y6, x6, w6);     
Psi(9) = direction (y9, x9, y6, x6, w6);     
Psi(10) = direction (y15, x15, y9, x9, w9);     
Psi(11) = direction (y1, x1, y9, x9, w9);     
Psi(12) = direction (y6, x6, y9, x9, w9);     
Psi(13) = direction (y1, x1, y15, x15, w15);     
Psi(14) = direction (y9, x9, y15, x15, w15); 
    
    check2=max(abs(L_hat-Psi'));
    
end

for i=4:length(X_hat)
    if X_hat(i)<0
        X_hat(i)=X_hat(i)+2*pi;
    end 
end

Adjusted_Unknown = X_0
Residuals=v
Adjusted_Observation=L_hat
    %Empirical reference standard deviation
    s_0 = sqrt(v'*P*v/r);
    %VC matrix of adjusted unknowns
    S_XX_hat =s_0^2*Q_xx;
    %Standard deviation of the adjusted unknows
    Std_Dev_Adj_unkwn =diag(sqrt(S_XX_hat))
    %Cofactor matrix of adjusted observations
    Q_LL_hat = A*Q_xx*A';
    %VC matrix of adjusted observations
    S_LL_hat = s_0^2*Q_LL_hat ;
    %Standard deviation of the adjusted observations
    Std_Dev_Adj_Obs = diag(sqrt(S_LL_hat))
    %Cofactor matrix of the residuals
    Q_vv = Q_LL - Q_LL_hat;
    %VC matrix of residuals
    S_vv = s_0^2*Q_vv;
    %Standard deviation of the residuals
    Std_Dev_Residuals = diag(sqrt(S_vv))




    
