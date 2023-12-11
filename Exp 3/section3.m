close all;
clear;
clc;

%% 3.3.b
w0 = pi * 0.15;
n = [0:199; 200:399; 400:599];
A = [2, 4, 0.5];
x = cos(w0 * n) .* A';
xn = x(1:end);
n = n(1:end);

plot(n, x(:)');

lambda = 0.9;
c0 = 0.5;
rho = 0.0;
b = 1- lambda;
a = [1 - lambda];

cn = filter(b, a, abs(xn));
gn = ones(1, length(cn));
gn(cn >= c0) = (cn(cn >= c0) / c0) .* (rho - 1);

yn = xn.*gn;

plot(xn);

% plot(n, cn);










