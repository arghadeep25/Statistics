%--------------------------------------------------------------------------
%   
%   SELECTED SECTIONS OF ADJUSTMENT CALCULATION 
%   t-Square Distribution
%   
%   Task           : 2        
%   Author         : Arghadeep Mazumder
%   Version        : May 23, 2017
%   
%--------------------------------------------------------------------------

clc;
clear all;
close all;

x=-5:0.0001:5;

A=tpdf(x,1);
figure;
plot(x,A,'r','DisplayName','t_1')
hold on

B=tpdf(x,2);
plot(x,B,'k','DisplayName','t_2')
hold on

C=tpdf(x,3);
plot(x,C,'b','DisplayName','t_3')
hold on

D=tpdf(x,4);
plot(x,D,'r-.','DisplayName','t_4')
hold on

E=tpdf(x,5);
plot(x,E,'k--','DisplayName','t_5')
hold on

F=tpdf(x,10);
plot(x,F,'b--','DisplayName','t_1_0')
hold on

G=tpdf(x,50);
plot(x,G,'k:','DisplayName','t_5_0')
hold on

H= normpdf(x,0,1);
plot(x,H,'m','DisplayName','N(0,1)')
hold on

title('t - Distribution')
xlabel('Values of Random Variable')
ylabel('Density Function')
legend('show')