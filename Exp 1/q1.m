close all;
clear;
clc;

% part 1.1
t = 0:0.01:2; % Define time axis
A = 5;
f = 1;
y = A*sin(2*pi*f*t); % Sine function

% Plotting the signals
figure('Name', 'Sine signal (Plot)');
plot(t, y, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;

figure('Name', 'Sine signal (Stem)');
stem(t, y, 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;


% part 1.2
t = 0:0.01:2; % Define time axis
A = 5;
f = 1;

y1 = A*sin(2*pi*f*t); % Normal sine

r = rand(1, length(y1)) - 0.5; % Noise
y2 = y1 + r; % Add noise to sine

% Plotting the signals
figure('Name', 'Normal Sine signal VS Noisy Sine Signal');
subplot(2,1,1);
plot(t,y1)
title('Normal Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(2,1,2);
plot(t,y2)
title('Noisy Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');

