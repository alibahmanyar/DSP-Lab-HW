%%
close all;
clear;
clc;

%% 2.1.a
% Testing myconv function against matlab's default conv
x = 1:2:10;
h = 1:5;
myconv(x,h)
conv(x,h)

%% 2.1.b
n = 0:200;
k = 50;
f = 1/k;

h = zeros(1, length(n)); % Integrator
h(1:10) = 0.1 * 1;

x = square(2 * pi * f * n) / 2 + 0.5;
y = myconv(h, x);

figure('Name', 'Integrator convolution');

subplot(2,1,1);
stem(n, x);
xlabel('n') ;
ylabel('x[n]') ;
title('Input Signal') ;
grid on ;

subplot(2,1,2);
stem(0:length(y) - 1, y);
xlabel('n') ;
ylabel('y[n]') ;
title('Output Signal') ;
grid on ;


%% 2.1.c
n = 0:14;
h2 = 0.25 * 0.75.^n;
y = myconv(h2, x);

figure('Name', 'Filter Convolution');
subplot(3,1,1);
stem(n, h2);
xlabel('n') ;
ylabel('h[n]') ;
title('Filter') ;
grid on ;

subplot(3,1,2);
stem(0:length(x) - 1, x);
xlabel('n') ;
ylabel('x[n]') ;
title('Input Signal') ;
xlim([0,200]);
grid on ;

subplot(3,1,3);
stem(0:length(y) - 1, y);
xlabel('n') ;
ylabel('y[n]') ;
title('Output Signal') ;
xlim([0,200]);
grid on ;


%% 2.1.d
% (0.2 (z^5 δ(n) - 5 z^4 δ(n) + 10 z^3 δ(n) - 10 z^2 δ(n) + 5 z δ(n) - δ(n)))/z^5
h3 = 0.2*[1,-5,10,-10,5,-1];
y = myconv(h3, x);

figure('Name', 'Filter Convolution');
subplot(3,1,1);
stem(h3);
xlabel('n') ;
ylabel('h[n]') ;
title('Filter') ;
grid on ;

subplot(3,1,2);
stem(0:length(x) - 1, x);
xlabel('n') ;
ylabel('x[n]') ;
title('Input Signal') ;
xlim([0,200]);
grid on ;

subplot(3,1,3);
stem(0:length(y) - 1, y);
xlabel('n') ;
ylabel('y[n]') ;
title('Output Signal') ;
xlim([0,200]);
grid on ;

%% 2.2.a
M  = 100;
n  = 0 : M - 1;
w1 = 0.05 * pi;
w2 = 0.20 * pi;
w3 = 0.35 * pi;
wa = 0.15 * pi;
wb = 0.25 * pi;

s  = sin(w2 * n);
v  = sin(w1 * n) + sin (w3 * n);
x  = s + v;

w  = 0.54 - 0.46 * sin(2 * pi * n / M);
h  = w .* ( wb * sinc(wb * (n - M/2)/pi)/pi - wa * sinc(wa * (n - M/2)/pi)/pi );

figure('Name', 'desired and final signal with noise');
subplot(2,1,1);
stem(n, x,'g');
title('X[n]'); xlabel('n'); ylabel(''); grid on;

subplot(2,1,2);
stem(n, s,'m');
title('S[n]'); xlabel('n'); ylabel(''); grid on;

%% 2.2.b
M  = 100;
n  = 0 : M - 1;
w1 = 0.05 * pi;
w2 = 0.20 * pi;
w3 = 0.35 * pi;
wa = 0.15 * pi;
wb = 0.25 * pi;

s  = sin(w2 * n);
v  = sin(w1 * n) + sin (w3 * n);
x  = s + v;

w  = 0.54 - 0.46 * sin(2 * pi * n / M);
h  = w .* ( wb * sinc(wb * (n - M/2)/pi)/pi - wa * sinc(wa * (n - M/2)/pi)/pi );

y = filter(h,1,x);

figure(5);
subplot(2,1,1);
stem(n, s,'r');
title('S[n]'); xlabel('n'); ylabel(''); grid on;

subplot(2,1,2);
stem(n, y,'k');
title('Y[n]'); xlabel('n'); ylabel(''); grid on;

%% 2.2.c
