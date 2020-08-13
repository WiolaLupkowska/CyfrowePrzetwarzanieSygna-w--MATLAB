close all
clear all

figure; 

for i=1:4
[RMS std_dev]= calculateRMS(0, 10, 0.1, 0, 10);
subplot(2,2,i); 
plot(std_dev,RMS);
xlabel('Standard Deviation');
ylabel('RMS');
grid minor;
hold on;
  
end
