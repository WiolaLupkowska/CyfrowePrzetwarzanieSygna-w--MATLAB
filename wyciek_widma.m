%close all
clear all;
clc
%some important parameters:
A   =   5;          %amplitude [a.u.]
f	=   3;          %Frquency [Hz]
phi =   0;          % Phase[rad]
ts  =   0;          %Starting time [s]
te  =   1.05;       %Ending time [s], bylo 1
fs  =   100;        %sampling rate [Hz], bylo 20
Ts  =   1/fs;       %sampling period [s]

t = 0 : Ts : te-Ts;             %Time domain [s], czemu ucinam o Ts??? transformata fouriera zanim policzy sygnal to replikuje, jeden sygnal zamkniety drugi otwarty na przedzialach lewo zamkniety, prawy otwarty
s = A * cos(2*pi*f * t + phi);  %signal [a.u.]
N = length(s);                  %signal length [samples]

w = taylorwin(N)'; %okno tyckeya, dawac rowne inne funkcje, rozne okienka,DO SPRAWKA!!!
%chebwin,gausswin,

%SPRAWKO
%wplyw funkcji okienkujacej na widmo
%wplyw cze prob na sygnal okienkowania
%wplyw okna na faze okienkowanego
%wplyw polozenia okna na widmo?

sw = s.*w;      %impose the window on the signal;

figure('Color','w','Position',[1e2 1e2 15e2 7e2])
subplot(311)
plot(t,s); hold on
plot(t, sw, 'r--'); hold off
legend('s','w')
xlabel('Time [s]'); ylabel('Signal [a.u.]'); grid;
xlim([t(1) t(end)])

%% Frequency domain analysis
S = fft(s) / N;
Sw = fft(sw) / N;
L= te - ts; %Signal length [s]
fline = (0 : 1/L : fs-1/L) - fs/2;  %dziedzina czestotliwosci 

%amplituda
subplot(312)
stem( fline, fftshift(abs(S))); hold on %duze S to fourier s, zeby spowrotem s to modul
stem(fline, fftshift(abs(Sw)), 'r--'); hold off
xlabel('Frequency [Hz]'); ylabel('Amplitude spectrum'); grid
xlim([min(fline) max(fline)])

%faza
subplot(313)
stem( fline, fftshift(angle(S))); hold on %duze S to fourier s, zeby spowrotem s to modul
stem(fline, fftshift(angle(Sw)), 'r--'); hold off
xlabel('Frequency [Hz]'); ylabel('Phase spectrum'); grid
xlim([min(fline) max(fline)])

