%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION 
%   t-Square datatribution
%   
%   Task           : 3        
%   Author         : Arghadeep Mazumder
%   Version        : May 23, 2017
%   
%--------------------------------------------------------------------------

clear;
clc;
close all;

data = load('distances.txt');
n = length(data); %Loading the file 'distances.txt'
dof = n-1; % Degree of Freedom
mv = mean(data); %Finding the Mean Value
Mean_Value= mv 
std_dev = std(data); %Calculating the Standard Deviation for Sigle Observaton
Standard_Deviation_for_Single_Obervation=std_dev
std_mean = std_dev/sqrt(n); %Calculating Standard Deviation for Arithmetic Mean
Standard_Deviation_for_Arithmetic_Mean=std_mean
res = ones(n,1)*mv - data;  %Calculating Residual Values
figure;
b=bar(res); % Plotting Residual Values in bar format
xlim([0 26]) % Setting the limits for the X-axis
title('Residual Values')
xlabel('Number of Observations')
ylabel('Residual Values [m]')

%Calculating Confidence Limit with S=95% for the expectation value as 
% well as for a single measurement

p_95 = .975;

a_95_so = mv - tinv(p_95,dof)*std_dev; 
b_95_so = mv + tinv(p_95,dof)*std_dev; 

Lower_Boundary_for_Single_Observation_95=a_95_so
Upper_Boundary_for_Single_Observation_95=b_95_so


a_95_mv = mv - tinv(p_95,dof)*std_mean; 
b_95_mv = mv + tinv(p_95,dof)*std_mean; 

Lower_Boundary_for_Mean_Value_95=a_95_mv
Upper_Boundary_for_Mean_Value_95=b_95_mv

%Calculating Confidence Limit with S=99% for the expectation value as 
% well as for a single measurement

p_99 = 0.995;
a_99_so = mv - tinv(p_99,dof)*std_dev; 
b_99_so = mv + tinv(p_99,dof)*std_dev; 

Lower_Boundary_for_Single_Observation_99=a_99_so
Upper_Boundary_for_Single_Observation_99=b_99_so


a_99_mv = mv - tinv(p_99,dof)*std_mean; 
b_99_mv = mv + tinv(p_99,dof)*std_mean; 

Lower_Boundary_for_Mean_Value_99=a_99_mv
Upper_Boundary_for_Mean_Value_99=b_99_mv
