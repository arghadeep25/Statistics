clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------
%Vector of Observations in rads
r31=206.9094*pi/200; 
r32=46.5027*pi/200;
r34=84.6449*pi/200; 
r35=115.5251*pi/200;
r36=155.5891*pi/200;
r=[r31;r32;r34;r35;r36];

L=[r34-r32; r35-r34; r36-r35; r31-r36];

%Error free_fixed Coordinates in m
y1=682.415; x1=321.052;
y2=203.526; x2=310.527;
y4=251.992; x4=506.222;
y5=420.028; x5=522.646;
y6=594.553; x6=501.494;
y=[y1 y2 y4 y5 y6]';
x=[x1 x2 x4 x5 x6]';

% Initial values for unknowns
x3=242;
y3=493;
X_0 = [x3; y3];

    %Number of observations
    no_r = length(r);
    no_n = length(L);

    %Number of unknowns
    no_u = length(X_0);

    %Redundancy
    r = no_n-no_u;
%--------------------------------------------------------------------------
%  stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sigma_r=0.001*pi/200;
S_LL = diag(ones(no_r,1)*sigma_r^2);

%VC matrix of the observations
%i=angles   j=directions
F=[0 -1 1 0 0;
   0 0 -1 1 0;
   0 0 0 -1 1;
   1 0 0 0 -1];

    S_LL=F*S_LL*F';

%Theoretical standard deviation
sigma_0 =1;

    %Cofactor matrix of the observations
    Q_LL =1/sigma_0^2*S_LL;

    %Weight matrix
    P =inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
    %break-off conditions
    epsilon = 10^-5;
    delta = 10^-12;
    max_x_hat = Inf;
    check2 = Inf;

    %Number of iterations
    iteration=0;

%Iteration
while max_x_hat>epsilon || check2>delta   % if sth is wrong and the iteration is endless...look at the design matrix for mistake

    %Vector of reduced observations
    t=[atan2((y1-y3),(x1-x3));
    atan2((y2-y3),(x2-x3));
    atan2((y4-y3),(x4-x3));
    atan2((y5-y3),(x5-x3));
    atan2((y6-y3),(x6-x3))];
    
L_0(1,1)=(t(3)-t(2));
L_0(2,1)=(t(4)-t(3));
L_0(3,1)=(t(5)-t(4));
L_0(4,1)=(t(1)-t(5));

l = L-L_0;

    %Designmatrix
dr_dx_342=dr_dx_angle(x3,y3,x4,y4,x2,y2);
dr_dx_354=dr_dx_angle(x3,y3,x5,y5,x4,y4);
dr_dx_365=dr_dx_angle(x3,y3,x6,y6,x5,y5);
dr_dx_316=dr_dx_angle(x3,y3,x1,y1,x6,y6);

dr_dy_342=dr_dy_angle(x3,y3,x4,y4,x2,y2);
dr_dy_354=dr_dy_angle(x3,y3,x5,y5,x4,y4);
dr_dy_365=dr_dy_angle(x3,y3,x6,y6,x5,y5);
dr_dy_316=dr_dy_angle(x3,y3,x1,y1,x6,y6);

A=[dr_dx_342 dr_dy_342;
    dr_dx_354 dr_dy_354;
    dr_dx_365 dr_dy_365;
    dr_dx_316 dr_dy_316];

    %Normal matrix
    N = A'*P*A;

    %Vector of absolute values
    n = A'*P*l;

    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);

    %Solution of normal equation
    x_hat = Q_xx*n;

    %Adjusted unknowns
    X_hat = X_0+x_hat;  

    %X_0 =>complete
    X_0 = X_hat;  

    %Update
x3 = X_0(1);
y3 = X_0(2);

    %Check 1
    % max_x_hat =>complete
    max_x_hat = max(abs(x_hat));

    %Vector of residuals
    v = A*x_hat-l;
    vTPv=v'*P*v;

    %Vector of adjusted observations
    L_hat = L+v;

    %Check 2
    %=>complete
    t=[atan2((y1-y3),(x1-x3));
    atan2((y2-y3),(x2-x3));
    atan2((y4-y3),(x4-x3));
    atan2((y5-y3),(x5-x3));
    atan2((y6-y3),(x6-x3))];

    Phi(1,1)=(t(3)-t(2));
    Phi(2,1)=(t(4)-t(3));
    Phi(3,1)=(t(5)-t(4));
    Phi(4,1)=(t(1)-t(5));
   
    check2 = max(abs(L_hat-Phi));

    %Update number of iterations
    iteration=iteration+1;

end
Adjusted_Observation=L_hat
Residuals=v
Adjusted_x3=x3
Adjusted_y3=y3

    %Empirical reference standard deviation
    s_0 = sqrt(v'*P*v/r);

    %VC matrix of adjusted unknowns
    S_XX_hat = s_0^2*Q_xx;

    %Standard deviation of the adjusted unknowns
    s_X = sqrt(diag(S_XX_hat));

    %Cofactor matrix of adjusted observations
    Q_LL_hat = A*Q_xx*A';

    %VC matrix of adjusted observations
    S_LL_hat = s_0^2*Q_LL_hat;

    %Standard deviation of the adjusted observations
    Std_dev_adjusted_unknown = sqrt(diag(S_LL_hat,0))

    %Cofactor matrix of the residuals
    Q_vv = Q_LL - Q_LL_hat;

    %VC matrix of residuals
    S_vv = s_0^2*Q_vv;

    %Standard deviation of the residuals
    Std_dev_Residuals = sqrt(diag(S_vv))
    
