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


%% 2.1.a
function y = myconv(h,x)
    M = length(h);
    L = length(x);
    N = M + L - 1;
    X = [x, zeros(1,M)];
    H = [h, zeros(1,L)];
    y = zeros(1, N);
    if(M > L)
        temp = X;
        X = H;
        H = temp;
    end
    for i = 1 : N
        for j = 1 : min(L, M)
            if (i - j + 1 > 0)
                y(i) = y(i) + H(j) * X(i - j + 1);
            end
        end
    end
end