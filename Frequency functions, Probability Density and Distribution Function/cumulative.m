% ##### Exercise 2 Adjustment Calculation #####

clc %clear all data in the command window
clear all %clear all data in the memory
close all %close all plots

data = load('distances.txt'); %Loading of the distance values from the memory
n=length(data); % Calculating the number of elements
sq=(sqrt(n)); % Calculating the number of Bins
[ybin,xbin]=hist(data,sq); 
bar(xbin,cumsum(ybin)/n*100) %Creating the histogram
hold on
plot(xbin,cumsum(ybin)/n*100,'r') %Plotting the polygon
xlim([min(data) max(data)]) % Limiting the X-axis
xlabel('Data','Fontsize',7) %To label the X-axis
ylabel('Cumulative Frequency','Fontsize',10) %To label the Y-axis