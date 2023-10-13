%%
close all;
clear;
clc;
%%
% part 1.1
t = 0:0.01:2;        % Define time axis
A = 5;               % Amplitude
f = 1;               % Frequency
y = A*sin(2*pi*f*t); % Sine function

% Plotting signals
figure('Name', 'Sine signal (Plot)');
plot(t, y, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;

figure('Name', 'Sine signal (Stem)');
stem(t, y,'color', '#7E2F8E', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Signal Amplitude');
title('Sine signal');
grid on;

%%
% part 1.2
r = rand(1, length(y)) - 0.5; % Noise
noised_signal = y + r;        % Add noise to sine, y defined in section 1.1

% Plotting signals
figure('Name', 'Normal Sine signal VS Noisy Sine Signal');
subplot(3,1,1);               % Plotting normal signal without noise
plot(t,y)
title('Normal Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(3,1,2);               % Plotting noisy signal
plot(t,noised_signal)
title('Noisy Signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
%hold on;
subplot(3,1,3);               % Plotting both signal together for comparing them
plot(t,y, t, noised_signal)
hold on;
title('Comparing signals');xlabel('Time (Second)');ylabel('Signal Amplitude');
legend('Normal signal', 'Noisy signal');
hold off;

%%
% part 1.3:
M1 = 0;
M2 = 20;
len = M2 + M1 + 1; % Window length

window = ones(1, len) / len; % Moving average window

y2 = conv(noised_signal, window);
t2 = 0:0.01:(0.01 * (length(y2) - 1));

plot(t2, y2);

% part 1.3.2:
M1 = 10;
M2 = 10;
len = M2 + M1 + 1; % Window length

window = ones(1, len) / len; % Moving average window

y3 = conv(noised_signal(M1:end), window);
t3 = 0:0.01:(0.01 * (length(y3) - 1));

plot(t3, y3);
legend("Noisy Signal", "Moving Avg M1=0, M2=20", "Moving Avg M1=10, M2=10");

%%
% part 1.4

