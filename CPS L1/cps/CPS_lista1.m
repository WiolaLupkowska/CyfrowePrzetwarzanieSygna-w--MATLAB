%% Zad.1
%T=1/f T=1 => f=1
t=-2:0.01:2;
y=15*cos(2*pi*1.*t+0);

figure;
plot(t,y);
grid minor;

%% List 1, Task 2
T = 0.01;  %Period
w = 2*pi/T;  %Frequency
A = 5;  %Amplitude

fi = pi/2; %Phase shift
tm = -fi/w;

t = [-3e-2 : e-4 : 3e-4];
xt = A*cos(w*t-fi);
figure('Color', 'w');
plot(t,xt);
%% Zad.3
w = 5*pi;
f = 2*pi/w;
fs = 100;  %Sampling frequency
t = [-4e-1 : 1/fs : 8e-1];

xt = 120*cos(5*pi.*t-pi/4); %Signal

figure;
plot(t,xt);
xlim([-0.4,0.8]);

xlabel('Amplitude [a.u.]')
ylabel('Time [s]')
%% Zad.4

t = 0 : 0.5 : 20;
xt = 7*cos(0.2*pi*t+0.5*pi);
figure;
plot(t,xt);
%T=2*pi/(0.2*pi);

%% Zad.5
t=-10:0.1:20;
w=pi/5;
xt=3*cos(w.*t-pi);

%figure
%plot(t,xt)

T=2*pi/w;

%%
caxis skaluje kolorki


%% Zad.6
w = 20*pi*f;
f = 10;
fs = 200;

t = 0.2 : 1/fs : 0.8;
x = 100 * cos(20 * pi * (t - 0.05));
% x=100cos(20pi*t-pi), z czego wynika ze fi=-pi

length(t) 
% ans=601, w tylu rownoodleg³ych punktach jest okreœlony sygna³
plot(t,x)
grid minor;
xlabel('Amplitude [a.u.]')
ylabel('Time [s]')
xlim([min(t),max(t)]);
%% Zad.7a
t = -10 : 0.1 : 10; T = 10;  A = 3;  f = 1/T;  t1 = -2;
tm = -2  %z wykresu
fi = -2*pi*f*tm  %=2*pi/5
fi_fake = pi/5;

x = A*cos(2*pi*f*t+fi_fake);
x1 = A*cos(2*pi*f*(t-t1));

figure
plot(t,x1,'-b',t,x,'or');


%NIE
%% Zad.7b
t = -10 : 0.1 : 10;  T = 10;  A = 3;  f = 1/T;
t1=5; tm=5; fi=-2*pi*f*tm;
fi_fake=pi;

x=A*cos(2*pi*f*t+fi_fake);
x1=A*cos(2*pi*f*(t-t1));

plot(t,x1)
hold on
plot(t,x)

%TAK
%% Zad.7c
t=-10:0.1:10; T=10; A=3; f=1/T;
t1=8; tm=-2; fi=-2*pi*f*tm;
fi_fake=2*pi/5;

x=A*cos(2*pi*f*t+fi_fake);
x1=A*cos(2*pi*f*(t-t1));

plot(t,x1)
hold on
plot(t,x)

%TAK

%% Zad.8

dt = 1/100;
t = -1:dt:1;
f0 = 2;
x = 300*real(exp(j*(2*pi*f0*(t-0.75))));
plot(t,x)
xlabel('Czas [s]'); ylabel('Amplituda [a.u.]'); 

x1=300*cos(4*pi*t+3*pi);
figure
plot(t,x1)
xlabel('Czas [s]'); ylabel('Amplituda [a.u.]'); 
%takie same

%% Zad.9
t=-5:0.1:10;
x1=4*exp(-j*pi/4);
x2=2*exp(-j*pi/2);
x3=x1+x2;
abs(x3);
angle(x3);

x3t=5.5959*cos(0.4*pi*t-1.0409);
%figure
%plot(t,x3t)

% z 7a
x=3*cos(pi/5*t+0.4*pi);
figure
plot(t,x,t,x3t)
% 

%na tym z 7 mieszcz¹ siê 2 pe³ne okresy
%na wykresie z tego zadania jeden

%% Zad.10

t=-10:0.1:10;
x1=5*exp(j*pi/3);
x2=7*exp(j*(-5*pi/4));
x3=3*exp(j*(3*pi/2));
x4=x1+x2+x3;

abs(x4)
angle(x4)
plot(t,x4)
%zle, na kartce ok
%wykres fazowy?

%% Zad.11
t=linspace(-10,20,300);
st=6*exp(-j*pi/3)*exp(j*pi/4*t);
s=imag(st);
%sit=6*j*sin(pi/4*t-pi/3);
sit=6*cos(pi/4*t+3*pi/3);


figure
plot(t,sit,'k')% Im
hold on
plot(t,st) %surowy sygnal
hold off
q=imag(diff(st)); %nie dzia³a bo sie length nie zgadza, ciagle length q sie zmiejsza o 1 jak sie zmieni w t
length(q)
figure
plot(t,q)


%% Zad.12
t=-10:0.1:10; T=10; A=3; f=1/T;
t1=5; tm=5; fi=-2*pi*f*tm;

x=A*cos(2*pi*f*t+fi);
%x1=A*cos(2*pi*f*(t-t1));
for i=1:length(t)
    k(i)=A*cos(2*pi*f*t(i)+fi);
end
    k;


wartSr=mean(x);
wartMax=max(x);
wartMin=min(x);
odchStd=std(x);

k(1)
k(2)
energia=0;
for i=1:length(t)
    energia=energia+k(i).^2;
    energia
    k(i)
end
energia

mocSr=mean(energia);
skuteczna=sqrt(mocSr);

%% Zad.13
n=1000;
y1=randn(n,1);
y2=rand(n,1);
t=linspace(0,10,n);
plot(t,y1,'rs')
hold on;
plot(t, y2,'ko')
hold off;
figure
hist(y1)
figure
hist(y2)


%% Zad.14
%lista 4 zadanie 5
t=linspace(0,1,1000);
sz=randn(1,length(t)); %œrednia jest 0 w normalnym)

%RMS=A/(sqrt(2));
RMS=[];
f_od_std=[];

for i=1:20
    f_od_std=[f_od_std i ];
    sz=sz*i;
    A= 0.5*(max(sz)-min(sz));
    RMS(i)=A/(sqrt(2));
end

figure;
plot(f_od_std,RMS)
hold off;

%% Zad.15
t=linspace(0,100,1000);
s=randn(1,length(t)).*8;  %moc szumu 1
%a=A*cos(80*pi*t);

A=0.5*(max(sygnal)-min(sygnal));%to jest odp
E=sum(s.^2);
Pn=mean(E);
PnDwa= (A^2)/2;
Ps=Pn*100; %to tez odp
A=sqrt(2*Ps);


%% Zad.16
t=linspace(0,20,100);
y=cos(t);
s=randn(1,length(t));
sygnal=y+s;

plot(t,y); %wykres cosinusa

En=sum(s.^2);
Pn=mean(En);
skut_sz=sqrt(Pn);
sr_n=mean(s);
abs_sr_n=abs(sr_n);
miedzyszczyt_n=max(s)-min(s);

Es=sum(y.^2);
Ps=mean(Es);
skut_sygnal=sqrt(Ps);
sr_s=mean(y);
sr_s=abs(sr_s);
miedzyszczyt_s=max(y)-min(y);


SNR=[];punkty=[];
for i=1:10
  punkty=[punkty 1*100]
  s=(randn(1,length(t)))*10*i^2;
  
  En=sum(s.^2);
  Pn=mean(En);
  snr=10*log10(Ps/Pn);
  SNR=[SNR snr];
  
end
SNR
figure
loglog(punkty,SNR);

%% Zad.17
t=-10:0.1:10;
x1=5*exp(j*pi/3);
x2=7*exp(j*(-5*pi/4));
x3=3*exp(j*(3*pi/2));
x4=x1+x2+x3;

abs(x4)
angle(x4)
plot(t,x4)


