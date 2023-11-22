close all;
clear;
clc;

%% 
N = 10e4;
f = [4, 8 ,12];
n = 0:N;

fs = 400;

t1 = 0:(1/fs):(2-1/fs);
t2 = 2:(1/fs):(4-1/fs);
t3 = 4:(1/fs):(6-1/fs);
x = [cos(2*pi*f(1)*t1),cos(2*pi*f(2)*t2),cos(2*pi*f(3)*t3)];



b1 = [0.969, -1.92377, 0.969531];
a1 = [1, -1.923772, 0.9399063];

b2 = [0.996088, -1.976468, 0.996088];
a2 = [1, -1.976468, 0.992177];

x_step = ones(1, N + 1);

y_step = filter(b1, a1, x_step);

figure('Name', 'Step Response');
stem(n, y_step, "LineWidth", 1.5);
title("Step Response of Filter H1 in Time Domain");
xlabel("n");
ylabel("Amplitude");
xlim([0 1000]);
grid on;

settling_time_value = find(abs(y_step - y_step(end)) >= 0.01, 1, 'last') + 1;

disp(sprintf("Setteling sample: %d", settling_time_value));

y_step = filter(b2, a2, x_step);

figure('Name', 'Step Response');
stem(n, y_step, "LineWidth", 1.5);
title("Step Response of Filter H2 in Time Domain");
xlabel("n");
ylabel("Amplitude");
xlim([0 1000]);
grid on;


settling_time_value = find(abs(y_step - y_step(end)) >= 0.01, 1, 'last') + 1;

disp(sprintf("Setteling sample: %d", settling_time_value));