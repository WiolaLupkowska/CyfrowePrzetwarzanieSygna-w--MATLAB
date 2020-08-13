
clear all;
clc

% t1,2      --- time
% N1,2      --- signal length
% s1,2      --- signal
% te1,2     --- ending time [s]
% ts        --- starting time [s]
% fs        --- sampling rate [Hz]
% Ts        --- sampling period [s]

ts = 0;
te1  =   1;       
te2  =   1.15; 
fs  =   50;      

Ts  =   1/fs; 
t1 = 0 : Ts : te1-Ts;            %odjêcie Ts: transformata Frourieta replikuje sygna³, przez co sygna³ musi byæ na jednym przedziale otwarty, na drugim zamkniêty 
t2 = 0 : Ts : te2-Ts; 
s1 = 2 * cos(2*pi*3 * t1 + pi);
s2 = 2 * cos(2*pi*3 * t2 + pi);
N1 = length(s1); 
N2 = length(s2);

w1 = hann(N1)';
w2 = hann(N2)';
sw1 = s1.*w1;   
sw2 = s2.*w2;   

figure('Color','w','Position',[1e2 1e2 15e2 5e2])
subplot(221)
plot(t1,s1); hold on
plot(t1, sw1, 'r--'); hold off
legend('s1','w1')
xlabel('Time1 [s]'); ylabel('Signal [a.u.]'); grid;
xlim([t1(1) t1(end)])

subplot(222)
plot(t2,s2); hold on
plot(t2, sw2, 'r--'); hold off
legend('s2','w2')
xlabel('Time2 [s]'); ylabel('Signal [a.u.]'); grid;
xlim([t2(1) t2(end)])

S1 = fft(s1) / N1;
Sw1 = fft(sw1) / N1
L1= te1-ts;                               %Signal length [s]
fline1 = (0 : 1/L1 : fs-1/L1) - fs/2;      %dziedzina czestotliwosci 

S2 = fft(s2) / N2;
Sw2 = fft(sw2) / N2;
L2= te2-ts;                               
fline2 = (0 : 1/L2 : fs-1/L2) - fs/2;      

%amplitude
subplot(223)
stem( fline1, fftshift(abs(S1))); hold on 
stem(fline1, fftshift(abs(Sw1)), 'r--'); hold off
xlabel('Frequency1 [Hz]'); ylabel('Amplitude spectrum'); grid
xlim([min(fline1) max(fline1)])

subplot(224)
stem( fline2, fftshift(abs(S2))); hold on 
stem(fline2, fftshift(abs(Sw2)), 'r--'); hold off
xlabel('Frequency2 [Hz]'); ylabel('Amplitude spectrum'); grid
xlim([min(fline2) max(fline2)])


