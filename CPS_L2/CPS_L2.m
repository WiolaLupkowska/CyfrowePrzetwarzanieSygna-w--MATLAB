%% list2

%% Task 1
t = 0:0.1:10;
t_norm = t./(max(t));
x = 8+10*cos(100*pi*t_norm+pi/3)+4*cos(200*pi*t_norm-pi/4);
x1 = 8+1/2*10*exp(1i*pi/3)*exp(1i*2*pi*50*t)+1/2*4*exp(-1i*pi/4)*exp(1i*2*pi*100*t)
+1/2*10*exp(-1i*pi/3)*exp(1i*2*pi*50*t)+1/2*4*exp(1i*pi/4)*exp(-1i*2*pi*100*t);

figure;
plot(t,x);
%% Task 2
%a)
t = 0:0.1:10;
t = t./max(t);
x = 10+4*exp(-1i*pi/4)*exp(1i*2*pi*10*t)+4*exp(1i*pi/4)*exp(-1i*2*pi*10*t)+2*exp(1i*pi/7)*exp(1i*2*pi*25*t)+2*exp(-1i*pi/7)*exp(-1i*2*pi*25*t);

figure;
plot(t,x);
%b)w=2*pi*k/T0;
f0 = 5 %bo f1=10, f2=25  => 2*pi*f0=w=2*pi/T =>
T0 = 1/f0

% c ) Ujemne czêstotliwoœci s¹ obecne w widmie, aby by³o ono symetryczne ~~Fourier?

%% Task 3
%a
t = 0:0.1:10;
t = t./max(t);
x = 20*cos(400*pi*t-pi/4)+5*cos(800*pi*t)-6*cos(1200*pi*t);
figure;
plot(t,x);
title("x(t)");
%b nie jest periodyczny
 
%c
y = x + 10 * cos( 600*pi*t + pi/6 );
figure;
plot(t,y);
title("y(t)");
%nie jest periodyczny
% ró¿ni¹ siê tym, ¿e maj¹ trochê inne szumy i nowy sygna³ ma amplitudê
% wiêksza o ok.10
 
%d
z = x + 10 *cos( 1800*pi*t + pi/6 );
figure;
plot(t,y);
title("z(t)");

%z na wykresie wszystkich sygnalow nachodzi na y
%zmienily sie szumy, amplitudy s¹ bardzo zbli¿one
figure;
plot(t,x,t,y,'+',t,z,'bd');
title("x(t), y(t), z(t)");

 %% Task 4
 t = 0:0.1:10;
t = t./max(t);
 
 w1 = 40*pi;
 w2 = 60*pi;
 w3 = 120*pi; %w=2*pi*f
 
  x= 2+4*cos(w1*t-pi/5)+3*sin(w2*t)+4*cos(w3*t-pi/3);
  
 f1 = 20;
 f2 = 30;
 f3 = 60; 
 
 f0 = 10;
 T0= 1/f0;
 
 w0=2*pi*f0;
 %x0=2+4*cos(40*pi*t-pi/5)+3*sin(60*pi*t)+4*cos(120*t-pi/3);
 xb = 2+4*exp(-1i*pi/5)*exp(1i*40*pi*t)+3*exp(1i*pi/2)*exp(1i*60*pi/3)+4*exp(-1i*pi/3)*exp(1i*120*pi*t);
 
%  b) widmo na krtce
% c)
w4=50*pi;
y = x + 10*cos(w4*t-pi/6);

figure;
plot(t,y,t,x);

% d) tak, jest periodyczny
w4=50*pi;
f4=25;
f0=5;
T0=1/f0;

%% Task 5

%% TASK 6
t=0:0.01:2;
fun=cos(21*pi*t).*sin(3*pi*t);
%fun= cos*(21*pi*t)*cos(3*pi*t-pi/2);
%fun=cos(21*pi*t+3*pit-pi/2)+cos(21*pi*t-3*pi*t+pi/2);
%fun=1/2*cos(24*pi*t-pi/2)+1/2*cos(18*pi*t+pi/2);
fun1=@(t)1/2*cos(24*pi*t-pi/2)*exp(-1i*(2*pi/T01)*k*t);
fun2=@(t)1/2*cos(18*pi*t+pi/2)*exp(-1i*(2*pi/T02)*k*t);
T01=1/12;
T02=1/9;

%plot(t,fun);
 k1= 0:1:T01;
 k2= 0:1:T02;
 
 
 
ak1 = 1/T01.*integral(fun1,0,12)
ak2 = 1/T02.*integral(fun1,0,12)





 
 