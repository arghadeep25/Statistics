clc
clear all
close all
data = load('distances.txt');
n=length(data);
bins=(sqrt(n));
hist(data,bins);
[yout,xout]=hist(data,bins);
%hold on
plot(xout,yout,'r')
bar(xout,cumsum(yout)/n*100)