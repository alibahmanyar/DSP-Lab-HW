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
r = rand(1, length(y)) - 0.5; % Noise
noised_signal = y + r; % Add noise to sine

% Plotting the signals
figure('Name', 'Normal Sine signal VS Noisy Sine Signal');
subplot(2,1,1);
plot(t,y)
title('Normal Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(2,1,2);
plot(t,noised_signal)
title('Noisy Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
hold on;

% part 1.3:
M1 = 0;
M2 = 20;
len = M2 - M1 + 1; % Window length

window = ones(1, len) / len; % Moving average window

y2 = conv(noised_signal, window);
t2 = 0:0.01:(0.01 * (length(y2) - 1));

plot(t2, y2);



