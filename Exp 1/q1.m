close all;
clear;
clc;

% part 1.1
t = 0:0.01:2; % Define time axis
A = 5;
f = 1;
y = A*sin(2*pi*f*t); % Sine function

% Plotting the signals
figure('Name', 'Sine signal');
plot(t, y, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;

figure('Name', 'Sine signal');
stem(t, y, 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;


