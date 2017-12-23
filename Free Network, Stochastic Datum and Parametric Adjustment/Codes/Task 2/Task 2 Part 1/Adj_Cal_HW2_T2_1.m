%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION
%                  HOMEWORK 2
%         Combined Horizontal Network 
% 
%   Author         : Arghadeep Mazumder
%   Mat. No        : 378554
%   Version        : May 30, 2017
%--------------------------------------------------------------------------


clc;
clear all;
close all;
format long g;

distances = load('Distances_1.txt');
directions = load('Directions_1.txt');
control_point = load('Control_Points.txt');
new_point = load('New_Points.txt');

%Gauss-Krueger coordinates for control points [m]
y1000 = control_point(1,2);
x1000 = control_point(1,3);
y2000 = control_point(2,2);
x2000 = control_point(2,3);
y3000 = control_point(3,2);
x3000 = control_point(3,3);
%New points [m]
y100 = new_point(1,2);
x100 = new_point(1,3);
y101 = new_point(2,2);
x101 = new_point(2,3);
y102 = new_point(3,2);
x102 = new_point(3,3);
y103 = new_point(4,2);
x103 = new_point(4,3);

%Initial values for orientation unknowns
w1000 = 0;
w2000 = 0;
w3000 = 0;
w100 = 0;
w101 = 0;
w102 = 0;
w103 = 0;


L = [distances(:,3); directions(:,3)*pi/200; y1000; x1000; y2000; x2000; y3000; x3000];    %Convert to [rad]


%Initial values for unknowns
X_0 = [y1000 x1000 y2000 x2000 y3000 x3000 y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';


%Number of observations
no_n = length(L);

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u; 

xy = reshape(X_0(1:6),1,6);

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
s_dist = 0.001;                       %[m]
s_dir = 0.001*pi/200;                 %Convert to [rad]
s_xy = 0.01;                          %[m]

s_LL = [s_dist^2*ones(length(distances),1); s_dir^2*ones(length(directions),1); s_xy^2*ones(length(xy),1)];
S_LL = diag(s_LL);

%Theoretical standard deviation
sigma_0 = 1;

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-7; %given accuracy in 0.001 
delta = 10^-9;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

%Initialising A
A = zeros(no_n,no_u);

%Iteration
while max_x_hat > epsilon

    %Vector of reduced distances
    L_0(1)=dis(y1000,x1000,y100,x100);
    L_0(2)=dis(y1000,x1000,y102,x102);
    L_0(3)=dis(y2000,x2000,y101,x101);
    L_0(4)=dis(y3000,x3000,y100,x100);
    L_0(5)=dis(y3000,x3000,y103,x103);
    L_0(6)=dis(y100,x100,y1000,x1000);
    L_0(7)=dis(y100,x100,y102,x102);
    L_0(8)=dis(y100,x100,y2000,x2000);
    L_0(9)=dis(y100,x100,y3000,x3000);
    L_0(10)=dis(y101,x101,y103,x103);
    L_0(11)=dis(y101,x101,y2000,x2000);
    L_0(12)=dis(y102,x102,y1000,x1000);
    L_0(13)=dis(y102,x102,y100,x100);
    L_0(14)=dis(y102,x102,y2000,x2000);
    L_0(15)=dis(y103,x103,y3000,x3000);
    L_0(16)=dis(y103,x103,y2000,x2000);
    L_0(17)=dis(y103,x103,y101,x101);
    
    %Vector of reduced directions
    L_0(18)=direction(y1000,x1000,y100,x100,w1000);
    L_0(19)=direction(y1000,x1000,y102,x102,w1000);
    L_0(20)=direction(y2000,x2000,y103,x103,w2000);
    L_0(21)=direction(y2000,x2000,y101,x101,w2000);
    L_0(22)=direction(y3000,x3000,y100,x100,w3000);
    L_0(23)=direction(y3000,x3000,y103,x103,w3000);
    L_0(24)=direction(y100,x100,y1000,x1000,w100);
    L_0(25)=direction(y100,x100,y2000,x2000,w100);
    L_0(26)=direction(y100,x100,y3000,x3000,w100);
    L_0(27)=direction(y101,x101,y103,x103,w101);
    L_0(28)=direction(y101,x101,y2000,x2000,w101);
    L_0(29)=direction(y102,x102,y1000,x1000,w102);
    L_0(30)=direction(y102,x102,y100,x100,w102);
    L_0(31)=direction(y102,x102,y2000,x2000,w102);
    L_0(32)=direction(y103,x103,y3000,x3000,w103);
    L_0(33)=direction(y103,x103,y2000,x2000,w103);
    L_0(34)=direction(y103,x103,y101,x101,w103);
    
    %Observed unknowns
    L_0(35) = y1000;
    L_0(36) = x1000;
    L_0(37) = y2000;
    L_0(38) = x2000;
    L_0(39) = y3000; 
    L_0(40) = x3000;
    
    
    l = L-L_0';

    %Design matrix
    A(1,1)=ds_dy_from(y1000,x1000,y100,x100); 
    A(1,2)=ds_dx_from(y1000,x1000,y100,x100); 
    A(1,7)=ds_dy_to(y1000,x1000,y100,x100);
    A(1,8)=ds_dx_to(y1000,x1000,y100,x100);
    
    A(2,1)=ds_dy_from(y1000,x1000,y102,x102); 
    A(2,2)=ds_dx_from(y1000,x1000,y102,x102);
    A(2,11)=ds_dy_to(y1000,x1000,y102,x102);
    A(2,12)=ds_dx_to(y1000,x1000,y102,x102);
    
    
    A(3,3)=ds_dy_from(y2000,x2000,y101,x101); 
    A(3,4)=ds_dx_from(y2000,x2000,y101,x101);
    A(3,9)=ds_dy_to(y2000,x2000,y101,x101);
    A(3,10)=ds_dx_to(y2000,x2000,y101,x101);

    A(4,5)=ds_dy_from(y3000,x3000,y100,x100);  
    A(4,6)=ds_dx_from(y3000,x3000,y100,x100);
    A(4,7)=ds_dy_to(y3000,x3000,y100,x100);
    A(4,8)=ds_dx_to(y3000,x3000,y100,x100);
    
    A(5,5)=ds_dy_from(y3000,x3000,y103,x103);  
    A(5,6)=ds_dx_from(y3000,x3000,y103,x103);
    A(5,13)=ds_dy_to(y3000,x3000,y103,x103);
    A(5,14)=ds_dx_to(y3000,x3000,y103,x103);
    
    A(6,7)=ds_dy_from(y100,x100,y1000,x1000);  
    A(6,8)=ds_dx_from(y100,x100,y1000,x1000);
    A(6,1)=ds_dy_to(y100,x100,y1000,x1000);
    A(6,2)=ds_dx_to(y100,x100,y1000,x1000);
    
    A(7,7)=ds_dy_from(y100,x100,y102,x102);  
    A(7,8)=ds_dx_from(y100,x100,y102,x102);
    A(7,11)=ds_dy_to(y100,x100,y102,x102);
    A(7,12)=ds_dx_to(y100,x100,y102,x102);
    
    A(8,7)=ds_dy_from(y100,x100,y2000,x2000);  
    A(8,8)=ds_dx_from(y100,x100,y2000,x2000);
    A(8,3)=ds_dy_to(y100,x100,y2000,x2000);
    A(8,4)=ds_dx_to(y100,x100,y2000,x2000);
    
    A(9,7)=ds_dy_from(y100,x100,y3000,x3000);  
    A(9,8)=ds_dx_from(y100,x100,y3000,x3000);
    A(9,5)=ds_dy_to(y100,x100,y3000,x3000);
    A(9,6)=ds_dx_to(y100,x100,y3000,x3000);
    
    A(10,9)=ds_dy_from(y101,x101,y103,x103);  
    A(10,10)=ds_dx_from(y101,x101,y103,x103);
    A(10,13)=ds_dy_to(y101,x101,y103,x103);
    A(10,14)=ds_dx_to(y101,x101,y103,x103);
    
    A(11,9)=ds_dy_from(y101,x101,y2000,x2000);  
    A(11,10)=ds_dx_from(y101,x101,y2000,x2000);
    A(11,3)=ds_dy_to(y101,x101,y2000,x2000);
    A(11,4)=ds_dx_to(y101,x101,y2000,x2000);
    
    A(12,11)=ds_dy_from(y102,x102,y1000,x1000);  
    A(12,12)=ds_dx_from(y102,x102,y1000,x1000);
    A(12,1)=ds_dy_to(y102,x102,y1000,x1000);
    A(12,2)=ds_dx_to(y102,x102,y1000,x1000);
    
    A(13,11)=ds_dy_from(y102,x102,y100,x100);  
    A(13,12)=ds_dx_from(y102,x102,y100,x100);
    A(13,7)=ds_dy_to(y102,x102,y100,x100);
    A(13,8)=ds_dx_to(y102,x102,y100,x100);
    
    A(14,11)=ds_dy_from(y102,x102,y2000,x2000);  
    A(14,12)=ds_dx_from(y102,x102,y2000,x2000);
    A(14,3)=ds_dy_to(y102,x102,y2000,x2000);
    A(14,4)=ds_dx_to(y102,x102,y2000,x2000);
    
    A(15,13)=ds_dy_from(y103,x103,y3000,x3000);  
    A(15,14)=ds_dx_from(y103,x103,y3000,x3000);
    A(15,5)=ds_dy_to(y103,x103,y3000,x3000);
    A(15,6)=ds_dx_to(y103,x103,y3000,x3000);
    
    A(16,13)=ds_dy_from(y103,x103,y2000,x2000);  
    A(16,14)=ds_dx_from(y103,x103,y2000,x2000);
    A(16,3)=ds_dy_to(y103,x103,y2000,x2000);
    A(16,4)=ds_dx_to(y103,x103,y2000,x2000);
    
    A(17,13)=ds_dy_from(y103,x103,y101,x101); 
    A(17,14)=ds_dx_from(y103,x103,y101,x101);
    A(17,9)=ds_dy_to(y103,x103,y101,x101);
    A(17,10)=ds_dx_to(y103,x103,y101,x101);	 
     
     
    A(18,1)=dr_dy_from(y1000,x1000,y100,x100); 
    A(18,2)=dr_dx_from(y1000,x1000,y100,x100); 
    A(18,7)=dr_dy_to(y1000,x1000,y100,x100);
    A(18,8)=dr_dx_to(y1000,x1000,y100,x100);
    A(18,15)=-1;
    
    A(19,1)=dr_dy_from(y1000,x1000,y102,x102); 
    A(19,2)=dr_dx_from(y1000,x1000,y102,x102);
    A(19,11)=dr_dy_to(y1000,x1000,y102,x102);
    A(19,12)=dr_dx_to(y1000,x1000,y102,x102);
    A(19,15)=-1;
    
    A(20,3)=dr_dy_from(y2000,x2000,y103,x103); 
    A(20,4)=dr_dx_from(y2000,x2000,y103,x103);
    A(20,13)=dr_dy_to(y2000,x2000,y103,x103);
    A(20,14)=dr_dx_to(y2000,x2000,y103,x103);
    A(20,16)=-1;
    
    A(21,3)=dr_dy_from(y2000,x2000,y101,x101); 
    A(21,4)=dr_dx_from(y2000,x2000,y101,x101);
    A(21,9)=dr_dy_to(y2000,x2000,y101,x101);
    A(21,10)=dr_dx_to(y2000,x2000,y101,x101);
    A(21,16)=-1;

    A(22,5)=dr_dy_from(y3000,x3000,y100,x100);  
    A(22,6)=dr_dx_from(y3000,x3000,y100,x100);
    A(22,7)=dr_dy_to(y3000,x3000,y100,x100);
    A(22,8)=dr_dx_to(y3000,x3000,y100,x100);
    A(22,17)=-1;
    
    A(23,5)=dr_dy_from(y3000,x3000,y103,x103);  
    A(23,6)=dr_dx_from(y3000,x3000,y103,x103);
    A(23,13)=dr_dy_to(y3000,x3000,y103,x103);
    A(23,14)=dr_dx_to(y3000,x3000,y103,x103);
    A(23,17)=-1;
    
    A(24,7)=dr_dy_from(y100,x100,y1000,x1000);  
    A(24,8)=dr_dx_from(y100,x100,y1000,x1000);
    A(24,1)=dr_dy_to(y100,x100,y1000,x1000);
    A(24,2)=dr_dx_to(y100,x100,y1000,x1000);
    A(24,18)=-1;
    
        
    A(25,7)=dr_dy_from(y100,x100,y2000,x2000);  
    A(25,8)=dr_dx_from(y100,x100,y2000,x2000);
    A(25,3)=dr_dy_to(y100,x100,y2000,x2000);
    A(25,4)=dr_dx_to(y100,x100,y2000,x2000);
    A(25,18)=-1;
    
    A(26,7)=dr_dy_from(y100,x100,y3000,x3000);  
    A(26,8)=dr_dx_from(y100,x100,y3000,x3000);
    A(26,5)=dr_dy_to(y100,x100,y3000,x3000);
    A(26,6)=dr_dx_to(y100,x100,y3000,x3000);
    A(26,18)=-1;
    
    A(27,9)=dr_dy_from(y101,x101,y103,x103);  
    A(27,10)=dr_dx_from(y101,x101,y103,x103);
    A(27,13)=dr_dy_to(y101,x101,y103,x103);
    A(27,14)=dr_dx_to(y101,x101,y103,x103);
    A(27,19)=-1;
    
    A(28,9)=dr_dy_from(y101,x101,y2000,x2000);  
    A(28,10)=dr_dx_from(y101,x101,y2000,x2000);
    A(28,3)=dr_dy_to(y101,x101,y2000,x2000);
    A(28,4)=dr_dx_to(y101,x101,y2000,x2000);
    A(28,19)=-1;
    
    A(29,11)=dr_dy_from(y102,x102,y1000,x1000);  
    A(29,12)=dr_dx_from(y102,x102,y1000,x1000);
    A(29,1)=dr_dy_to(y102,x102,y1000,x1000);
    A(29,2)=dr_dx_to(y102,x102,y1000,x1000);
    A(29,20)=-1;
    
    A(30,11)=dr_dy_from(y102,x102,y100,x100); 
    A(30,12)=dr_dx_from(y102,x102,y100,x100);
    A(30,7)=dr_dy_to(y102,x102,y100,x100);
    A(30,8)=dr_dx_to(y102,x102,y100,x100);
    A(30,20)=-1;
    
    A(31,11)=dr_dy_from(y102,x102,y2000,x2000);  
    A(31,12)=dr_dx_from(y102,x102,y2000,x2000);
    A(31,3)=dr_dy_to(y102,x102,y2000,x2000);
    A(31,4)=dr_dx_to(y102,x102,y2000,x2000);
    A(31,20)=-1;
    
    A(32,13)=dr_dy_from(y103,x103,y3000,x3000);  
    A(32,14)=dr_dx_from(y103,x103,y3000,x3000);
    A(32,5)=dr_dy_to(y103,x103,y3000,x3000);
    A(32,6)=dr_dx_to(y103,x103,y3000,x3000);
    A(32,21)=-1;
    
    A(33,13)=dr_dy_from(y103,x103,y2000,x2000);  
    A(33,14)=dr_dx_from(y103,x103,y2000,x2000);
    A(33,3)=dr_dy_to(y103,x103,y2000,x2000);
    A(33,4)=dr_dx_to(y103,x103,y2000,x2000);
    A(33,21)=-1;
    
    A(34,13)=dr_dy_from(y103,x103,y101,x101);  
    A(34,14)=dr_dx_from(y103,x103,y101,x101);
    A(34,9)=dr_dy_to(y103,x103,y101,x101);
    A(34,10)=dr_dx_to(y103,x103,y101,x101); 
    A(34,21)=-1;
    
    A(35,1)=1;
    A(36,2)=1;
    A(37,3)=1;
    A(38,4)=1;
    A(39,5)=1;
    A(40,6)=1;
    
    
    %Normal matrix
    N = A'*P*A;

        
    %Vector of right hand side of normal equations
     n = A'*P*l;
    
    %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);
    
    
    %Solution of normal equation
     x_hat = Q_xx*n;
    
    %Adjusted unknowns
     X_hat = X_0+x_hat;
    
    %Update
    X_0 = X_hat;
    
    y1000 = X_0(1);
    x1000 = X_0(2);
    y2000 = X_0(3);
    x2000 = X_0(4);
    y3000 = X_0(5);
    x3000 = X_0(6);
    y100 = X_0(7);
    x100 = X_0(8);
    y101 = X_0(9);
    x101 = X_0(10);
    y102 = X_0(11);
    x102 = X_0(12);
    y103 = X_0(13);
    x103 = X_0(14);
    w1000 = X_0(15);
    w2000 = X_0(16);
    w3000 = X_0(17);
    w100 = X_0(18);
    w101 = X_0(19);
    w102 = X_0(20);
    w103 = X_0(21);
    
    %Check 1
     max_x_hat = max(abs(x_hat));
    
    %Update number of iterations
    iteration=iteration+1;

end

%Convert to [gon] and check the quadrants
gon = X_0(15:21,1)*200/pi;     %convert w to gon

gon = gon+400;

%gon(1) = gon(1)+400;          %check w1-w9 in X_0 >> check the quadrants manually >> add +400 for the minus value
%gon(2) = gon(2)+400;
%gon(3) = gon(3)+400;
%gon(4) = gon(4)+400;




%Vector of residuals
v = A*x_hat-l;
v_gon = v(18:34,1)*200/pi;             %Convert to [gon]

%Objective function
vTPv = v'*P*v;

%Vector of adjusted observations
L_hat = L+v;
L_hat_gon = L_hat(18:34)*200/pi;      %Convert to [gon]


%Final check for the linearization
    
    %Vector of reduced distances
    Phi_X_hat(1)=dis(y1000,x1000,y100,x100);
    Phi_X_hat(2)=dis(y1000,x1000,y102,x102);
    Phi_X_hat(3)=dis(y2000,x2000,y101,x101);
    Phi_X_hat(4)=dis(y3000,x3000,y100,x100);
    Phi_X_hat(5)=dis(y3000,x3000,y103,x103);
    Phi_X_hat(6)=dis(y100,x100,y1000,x1000);
    Phi_X_hat(7)=dis(y100,x100,y102,x102);
    Phi_X_hat(8)=dis(y100,x100,y2000,x2000);
    Phi_X_hat(9)=dis(y100,x100,y3000,x3000);
    Phi_X_hat(10)=dis(y101,x101,y103,x103);
    Phi_X_hat(11)=dis(y101,x101,y2000,x2000);
    Phi_X_hat(12)=dis(y102,x102,y1000,x1000);
    Phi_X_hat(13)=dis(y102,x102,y100,x100);
    Phi_X_hat(14)=dis(y102,x102,y2000,x2000);
    Phi_X_hat(15)=dis(y103,x103,y3000,x3000);
    Phi_X_hat(16)=dis(y103,x103,y2000,x2000);
    Phi_X_hat(17)=dis(y103,x103,y101,x101);
    
    %Vector of reduced directions
    Phi_X_hat(18)=direction(y1000,x1000,y100,x100,w1000);
    Phi_X_hat(19)=direction(y1000,x1000,y102,x102,w1000);
    Phi_X_hat(20)=direction(y2000,x2000,y103,x103,w2000);
    Phi_X_hat(21)=direction(y2000,x2000,y101,x101,w2000);
    Phi_X_hat(22)=direction(y3000,x3000,y100,x100,w3000);
    Phi_X_hat(23)=direction(y3000,x3000,y103,x103,w3000);
    Phi_X_hat(24)=direction(y100,x100,y1000,x1000,w100);
    Phi_X_hat(25)=direction(y100,x100,y2000,x2000,w100);
    Phi_X_hat(26)=direction(y100,x100,y3000,x3000,w100);
    Phi_X_hat(27)=direction(y101,x101,y103,x103,w101);
    Phi_X_hat(28)=direction(y101,x101,y2000,x2000,w101);
    Phi_X_hat(29)=direction(y102,x102,y1000,x1000,w102);
    Phi_X_hat(30)=direction(y102,x102,y100,x100,w102);
    Phi_X_hat(31)=direction(y102,x102,y2000,x2000,w102);
    Phi_X_hat(32)=direction(y103,x103,y3000,x3000,w103);
    Phi_X_hat(33)=direction(y103,x103,y2000,x2000,w103);
    Phi_X_hat(34)=direction(y103,x103,y101,x101,w103);
    
    % Control points
    Phi_X_hat(35) = y1000;
    Phi_X_hat(36) = x1000;
    Phi_X_hat(37) = y2000;
    Phi_X_hat(38) = x2000;
    Phi_X_hat(39) = y3000; 
   Phi_X_hat(40) = x3000;

% Final Check
Finalcheck = max(abs(L_hat-Phi_X_hat'));

if Finalcheck<=delta
    disp('everything is fine!')
else
    disp('Something is wrong!')
end





%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
s_X_gon = s_X(15:21,1)*200/pi;        %Convert to [gon]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
s_L_hat_gon = s_L_hat(18:34)*200/pi;    %Convert to [gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
s_v_gon = s_v(18:34,1)*200/pi;        %Convert to [gon]



%--------------------------------------------------------------------------
%  Global test
%--------------------------------------------------------------------------

%% Global test
Tx2 = (r*s_0^2)/sigma_0^2;

%for S = 95%
tx2u = chi2inv(0.975,r);  %1-a/2
tx2l = chi2inv(0.025,r);  %a/2

if tx2l<Tx2 && Tx2<tx2u
  disp('Fails to reject the H0')
else
  disp('Rejects the H0')
end


%--------------------------------------------------------------------------
% Internal and external reliability parameters
%--------------------------------------------------------------------------
% Parameters for internal
% Redundancy numbers
EV = diag(Q_vv*P);  %how much error has been transfered to residuals
EV_100 = EV*100;

% Standardised residuals
sigma_v = sigma_0^2*sqrt(diag(Q_vv));  %theoritical std dev that has no blunders
NV = abs(v)./sigma_v;

% Potential magnitude of a blunder
GF = -v./(diag(Q_vv*P));
GF_gon = GF(18:34,1)*200/pi;  %only the obs. directions

% Lower boundary value for blunders
GRZW = ones(no_n,1);

for i = 1:no_n
  GRZW(i,1) = sigma_0*4.13/(sqrt(EV(i,1)*P(i,i)));
end
GRZW_d= GRZW(1:17,1);
GRZW_gon = GRZW(18:34,1)*200/pi;

%Parameters for external reliability
P_diag = diag(P);
r_w = ones(17,1);

r_w(1) = P_diag(18,1)/(P_diag(18,1)+P_diag(19,1));
r_w(2) = P_diag(19,1)/(P_diag(18,1)+P_diag(19,1));
r_w(3) = P_diag(20,1)/(P_diag(20,1)+P_diag(21,1));
r_w(4) = P_diag(21,1)/(P_diag(20,1)+P_diag(21,1));
r_w(5) = P_diag(22,1)/(P_diag(22,1)+P_diag(23,1));
r_w(6) = P_diag(23,1)/(P_diag(22,1)+P_diag(23,1));
r_w(7) = P_diag(24,1)/(P_diag(24,1)+P_diag(25,1)+P_diag(26,1));
r_w(8) = P_diag(25,1)/(P_diag(24,1)+P_diag(25,1)+P_diag(26,1));
r_w(9) = P_diag(26,1)/(P_diag(24,1)+P_diag(25,1)+P_diag(26,1));
r_w(10) = P_diag(27,1)/(P_diag(27,1)+P_diag(28,1));
r_w(11) = P_diag(28,1)/(P_diag(27,1)+P_diag(28,1));
r_w(12) = P_diag(29,1)/(P_diag(29,1)+P_diag(30,1)+P_diag(31,1));
r_w(13) = P_diag(30,1)/(P_diag(29,1)+P_diag(30,1)+P_diag(31,1));
r_w(14) = P_diag(31,1)/(P_diag(29,1)+P_diag(30,1)+P_diag(31,1));
r_w(15) = P_diag(32,1)/(P_diag(32,1)+P_diag(33,1)+P_diag(34,1));
r_w(16) = P_diag(33,1)/(P_diag(32,1)+P_diag(33,1)+P_diag(34,1));
r_w(17) = P_diag(34,1)/(P_diag(32,1)+P_diag(33,1)+P_diag(34,1));

%dd = dist(:,3);

dd(1) = distances(1,3);
dd(2) = distances(2,3);
dd(3) = distances(15,3);
dd(4) = distances(3,3);
dd(5) = distances(4,3);
dd(6) = distances(5,3);
dd(7) = distances(6,3);
dd(8) = distances(8,3);
dd(9) = distances(9,3);
dd(10) = distances(10,3);
dd(11) = distances(11,3);
dd(12) = distances(12,3);
dd(13) = distances(13,3);
dd(14) = distances(14,3);
dd(15) = distances(15,3);
dd(16) = distances(16,3);
dd(17) = distances(17,3);

%Impact of the boundary value on the coordinates of the corresponding
%points
EGK = ones(40,1); 
for i=1:17
  EGK(i) = (1-EV(i,1))*GRZW(i,1);
end

for i=18:34
  EGK(i) = (1-EV(i,1)-r_w(i-17,1))*GRZW(i,1)*dd(i-17);
end

for i=35:40 
  EGK(i) = (1-EV(i,1))*GRZW(i,1);
end
  

EGK_gon = EGK(18:34,1)*200/pi;

%Impact of a potential blunder on a point corresponding to the measurement
EP = ones(40,1); 

for i=1:17
  EP(i,1) = (1-EV(i,1))*GF(i,1);
end

for i=18:34
  EP(i,1) = (1-EV(i,1)-r_w(i-17,1))*GF(i,1)*dd(i-17);
end

for i=35:40 
  EP(i,1) = (1-EV(i,1))*GF(i,1);
end

EP_gon = EP(18:34,1)*200/pi;


results.L=[L(1:18,1); directions(:,3)];
results.v=[v(1:18,1); v_gon(:,1)];
results.L=[L_hat(1:18,1); L_hat_gon(:,1)];
results.s=[s_v(1:18,1); s_v_gon(:,1)];
results.s=[s_L_hat(1:18,1); s_L_hat_gon(:,1)];
results.GF=[GF(1:18,1); GF_gon(:,1)];
results.GRZW=[GRZW(1:18,1); GRZW_gon(:,1)];
results.EGK=[EGK(1:18,1); EGK_gon(:,1)];
results.EP=[EP(1:18,1); EP_gon(:,1)];

results= struct2table(results);
writetable(results, 'task2_1.xls');