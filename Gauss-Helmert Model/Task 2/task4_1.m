clc;
clear all;
close all;
format long g

% Break-off condition and initial settings
ep = 1e-12;
c = 1e-5;
max_delta_v = inf;
epsilon=10^-12;
delta= 10^-12;

%Number of iterations
it = 0;
iteration=0;

%--------------------------------------------------------------------------
%   Observations and initial values for unknowns
%--------------------------------------------------------------------------

data=load('plane.txt');

x=data(:,2);
y=data(:,3);
z=data(:,4);
L=[x;y;z];

%Initial values for unknowns
nx=1;
ny=1;
nz=1;
d=1;

n_v=[nx, ny, nz]';

X_0=[n_v;d];

%Initial values for the residuals
vx= zeros(length(x),1);
vy= zeros(length(y),1);
vz= zeros(length(z),1);

v=[vx;vy;vz];

%Number of condition equations
no_n=length(data);

%Number of unknowns
no_u=length(X_0);

%Redundancy
r=no_n-no_u+1;

%--------------------------------------------------------------------------
%  Initial stochastic model
%--------------------------------------------------------------------------
%Weight matrix
P= eye(3*no_n,3*no_n);
Q_ll=inv(P);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------

while max_delta_v>ep
 v_appr=v;
 max_PsiU=inf;
 max_X_hat=inf;
 while or(max_X_hat>delta,max_PsiU>epsilon)

    %Condition equations Psi_i
    Psi= nx.*(x+vx)+ny.*(y+vy)+nz.*(z+vz)-d;

    %Designmatrix A
    A=[x+vx y+vy z+vz -ones(no_n,1)];
    
    %Designmatrix B
    B1=nx*eye(no_n,no_n);
    B2=ny*eye(no_n,no_n);
    B3=nz*eye(no_n,no_n);
    B=[B1, B2, B3];
    %Designmatrix C
    C= [2*nx 2*ny 2*nz 0];
    
    %Vector of misclosures
    w_1= -B*v+Psi;
    w_2= nx^2+ny^2+nz^2-1;
    
    %Normal matrix
    N= [-A'*(B*Q_ll*B')^-1*A C'; C 0]; 

    %Vector of the right hand side
    n= [A'*(B*Q_ll*B')^-1*w_1; -w_2];

    %Inversion of normal matrix
    N_inv=inv(N);
    
    %Solution of normal equation
    x_hat=N_inv*n;
    
    %Adjusted unknowns
    X_hat= X_0+x_hat(1:end-1);
    
    %Update of the unknowns
    X_0=X_hat;
    
	nx=X_0(1);
    ny=X_0(2);
    nz=X_0(3);
    d=X_0(4);
    
    %Lagrangian Multipliers
    k_1= -(B*Q_ll*B')^-1*(A*x_hat(1:end-1)+w_1);
    
    %Residuals
    v= Q_ll*B'*k_1;
    
    %Update of the residuals
    vx= v(1:no_n);
    vy= v(no_n+1:2*no_n);
    vz= v(2*no_n+1:end);
    
    %Check 1
    max_X_hat = max(abs(x_hat));
    %Check 2
    PsiU = nx.*(x+vx)+ny.*(y+vy)+nz.*(z+vz)-d;
    max_PsiU = max(abs(PsiU));
    
    %Update number of iterations
    iteration=iteration+1;

end
    
    %Update of the stochastic model
    P=  diag(1./(abs(v)+c));
    Q_ll=inv(P);
    
    %Check
    max_delta_v= max(abs(v-v_appr));
    
    %Update number of iterations
    it=it+1;

end

%Vector of adjusted observations
L_hat=L+v;

%Empirical reference standard deviation
s_0=sqrt(v'*inv(Q_ll)*v/r);

%VC matrix of adjusted unknowns
Q_xx= -N_inv(1:no_u,1:no_u);
S_XX_hat=s_0^2*Q_xx;

%Standard deviation of the adjusted unknows
s_X=sqrt(diag(S_XX_hat));

%Cofactor matrix of the residuals
Q_vv= Q_ll*B'*(B*Q_ll*B')^-1*B*Q_ll;

%VC matrix of residuals
S_vv=s_0^2*Q_vv;

%Standard deviation of the residuals
s_v=sqrt(diag(S_vv));

%Cofactor matrix of adjusted observations
Q_LL_hat=Q_ll-Q_vv;

%VC matrix of adjusted observations
S_LL_hat=s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat=sqrt(diag(S_LL_hat));

figure
bar(v*1000)
xlabel('Number of Residual')
ylabel('v (mm)')
title('Residuals')
result.L=[L(1:30,1)];
% result.x=x;
result.s_L_hat=[s_L_hat(1:30,1)];
result.v=[v(1:30,1)];
result.s_v=[s_v(1:30,1)];
result= struct2table(result);
writetable(result, 'task2_1.xls');

