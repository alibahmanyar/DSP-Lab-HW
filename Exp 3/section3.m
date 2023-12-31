%%
close all;
clear;
clc;
%% 3.3.a
fs = 2000;
t = 0:1 / fs:1;
A = 10;
x = A * cos(2 * pi * 50 * t);

disp("MS of X:");
ms_val_x = mean(x .^ 2);
disp(ms_val_x);
disp(A ^ 2/2);
disp(rms(x) ^ 2);

disp('Mean of x');
disp(mae(x));
disp(2 * A / pi);

%% Input signal
w0 = pi * 0.15;
n = [0:199; 200:399; 400:599];
A = [2 4 0.5];
x = cos(w0 * n) .* A';
x = x';
x = x(1:end);
n = n';
n = n(1:end);
%% 3.3.b
% Control Signal
lambda = 0.9;
c0 = 0.5;
rho = 0.2;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

% Gain Signal
gn = ones(1, length(cn));
gn(cn >= c0) = (cn(cn >= c0) / c0) .^ (rho - 1);
yn = gn .* x;


% Plot Signals
% fig 3.6
figure("Name", 'Input Signal vs Compressed Signal (3.6)');
subplot(1, 2, 1);
plot(n, x, 'LineWidth', 1.5);
title('Input Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, yn, 'LineWidth', 1.5);
title('Compressed Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;

% fig 3.7
figure('Name', 'Contol Signal vs Gain Signal (3.7)');
subplot(1, 2, 1);
plot(n, cn, 'LineWidth', 1.5);
title('Control Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn, 'LineWidth', 1.5);
title('Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 1.2]);
grid on;

L = 7;
cn_smoothed = movmean(cn, [L - 1 0]);
gn_smoothed = ones(1, length(cn_smoothed));
gn_smoothed(cn_smoothed >= c0) = (cn_smoothed(cn_smoothed >= c0) / c0) .^ (rho - 1);
yn_smoothed = gn_smoothed .* x;

% fig 3.8
figure('Name', 'Smoothed Output Signal vs Smoothed Gain Signal (3.8)');
subplot(1, 2, 1);
plot(n, yn_smoothed, 'LineWidth', 1.5);
title('Smoothed Output Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn_smoothed, 'LineWidth', 1.5);
title('Smoothed Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-2 2]);
grid on;

% fig.3.9
lambda = 0.1;
c0 = 1.5;
rho = 0.1;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

L = 7;
cn_smoothed = movmean(cn, L);
gn_smoothed = ones(1, length(cn_smoothed));
gn_smoothed(cn_smoothed >= c0) = (cn_smoothed(cn_smoothed >= c0) / c0) .^ (rho - 1);
yn_smoothed = gn_smoothed .* x;

figure('Name', 'Smoothed Output Signal vs Smoothed Gain Signal (3.9)');
subplot(1, 2, 1);
plot(n, yn_smoothed, 'LineWidth', 1.5);
title('Smoothed Output Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn_smoothed, 'LineWidth', 1.5);
title('Smoothed Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-2 2]);
grid on;

%% fig 3.10
lambda = 0.9;
c0 = 0.5;
rho = 2;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

L = 7;
cn_smoothed = movmean(cn, L);
gn_smoothed = ones(1, length(cn_smoothed));
gn_smoothed(cn_smoothed <= c0) = (cn_smoothed(cn_smoothed <= c0) / c0) .^ (rho - 1);
yn_smoothed = gn_smoothed .* x;

figure('Name', 'Smoothed Output Signal vs Smoothed Gain Signal (3.10)');
subplot(1, 2, 1);
plot(n, yn_smoothed, 'LineWidth', 1.5);
title('Smoothed Output Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn_smoothed, 'LineWidth', 1.5);
title('Smoothed Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-2 2]);
grid on;

%% fig 3.11
lambda = 0.1;
c0 = 0.5;
rho = 10;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

L = 7;
cn_smoothed = movmean(cn, L);
gn_smoothed = ones(1, length(cn_smoothed));
gn_smoothed(cn_smoothed <= c0) = (cn_smoothed(cn_smoothed <= c0) / c0) .^ (rho - 1);
yn_smoothed = gn_smoothed .* x;

figure('Name', 'Smoothed Output Signal vs Smoothed Gain Signal (3.10)');
subplot(1, 2, 1);
plot(n, yn_smoothed, 'LineWidth', 1.5);
title('Smoothed Output Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn_smoothed, 'LineWidth', 1.5);
title('Smoothed Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-2 2]);
grid on;

%% 3.3.c
% Control Signal
lambda = 0.9;
c0 = 0.5;
rho = 0.25;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

% Gain Signal
gn = ones(1, length(cn));
gn(cn >= c0) = (cn(cn >= c0) / c0) .^ (rho - 1);
yn = gn .* x;


% Plot Signals
% fig 3.6
figure("Name", 'Input Signal vs Compressed Signal (3.6)');
subplot(1, 2, 1);
plot(n, x, 'LineWidth', 1.5);
title('Input Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, yn, 'LineWidth', 1.5);
title('Compressed Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;

% fig 3.7
figure('Name', 'Contol Signal vs Gain Signal (3.7)');
subplot(1, 2, 1);
plot(n, cn, 'LineWidth', 1.5);
title('Control Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn, 'LineWidth', 1.5);
title('Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 1.2]);
grid on;

%% 3.3.d
% Control Signal
lambda = 0.9;
c0 = 1.3;
rho = 0.5;
b = 1 - lambda;
a = [1, -lambda];
cn = filter(b, a, abs(x));

% Gain Signal
gn = ones(1, length(cn));
gn(cn >= c0) = (cn(cn >= c0) / c0) .^ (rho - 1);
yn = gn .* x;


% Plot Signals
% fig 3.6
figure("Name", 'Input Signal vs Compressed Signal (3.6)');
subplot(1, 2, 1);
plot(n, x, 'LineWidth', 1.5);
title('Input Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;
subplot(1, 2, 2);
plot(n, yn, 'LineWidth', 1.5);
title('Compressed Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([-5 5]);
grid on;

% fig 3.7
figure('Name', 'Contol Signal vs Gain Signal (3.7)');
subplot(1, 2, 1);
plot(n, cn, 'LineWidth', 1.5);
title('Control Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 5]);
grid on;
subplot(1, 2, 2);
plot(n, gn, 'LineWidth', 1.5);
title('Gain Signal');
xlabel('n');
ylabel('amplitude');
xlim([0 600]);
ylim([0 1.2]);
grid on;