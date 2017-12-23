%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION 
%   Chi-Square Distribution
% 
%   Task           : 2        
%   Author         : Arghadeep Mazumder
%   Version        : May 23, 2017
%   
%--------------------------------------------------------------------------

clc;
clear all;
close all;

x=0:0.20:15;

A=chi2pdf(x,1);
figure;
plot(x,A,'r','DisplayName','x^2_f_1')
hold on

B=chi2pdf(x,2);
plot(x,B,'k','DisplayName','x^2_f_2')
hold on

C=chi2pdf(x,3);
plot(x,C,'b','DisplayName','x^2_f_3')
hold on

D=chi2pdf(x,4);
plot(x,D,'r--','DisplayName','x^2_f_4')
hold on

E=chi2pdf(x,5);
plot(x,E,'k--','DisplayName','x^2_f_5')
hold on

F=chi2pdf(x,7);
plot(x,F,'b--','DisplayName','x^2_f_7')
hold on

G=chi2pdf(x,9);
plot(x,G,'k:','DisplayName','x^2_f_9')
hold on

title('X^2 - Distribution')
xlabel('Values of Random Variable')
ylabel('Density Function')
legend('show')
