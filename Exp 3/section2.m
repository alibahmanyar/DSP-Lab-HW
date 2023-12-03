close all;
clear;
clc;

%%
% H1:
b1 = [0.969531, -1.923772, 0.969531];
a1 = [1, -1.923772, 0.939063];

% H2:
b2 = [0.996088, -1.976468, 0.996088];
a2 = [1, -1.976468, 0.992177];
%% 3-2-a
% H1:
figure('name', 'H1 response')
freqz(b1,a1,1024,'whole');
% H2:
figure('name', 'H2 response')
freqz(b2,a2,1024,'whole');

% Canonical Presentation
csys1  = canon(tf(b1 , a1),'companion')
csys2  = canon(tf(b2 , a2),'companion')

%% 3.2.b
N = 10e4;
n = 0:N;
x_step = ones(1, N + 1);

% H1
y_step = filter(b1, a1, x_step);

figure('Name', 'Step Response');
plot(n, y_step, "LineWidth", 1.5);
title("Step Response of Filter H1 in Time Domain");
xlabel("n");
ylabel("Amplitude");
xlim([0 1000]);
grid on;

settling_time_value = find(abs(y_step - y_step(end)) >= 0.01, 1, 'last') + 1;
fprintf("Setteling sample (H1): %d\n", settling_time_value);

% H2
y_step = filter(b2, a2, x_step);

figure('Name', 'Step Response');
plot(n, y_step, "LineWidth", 1.5);
title("Step Response of Filter H2 in Time Domain");
xlabel("n");
ylabel("Amplitude");
xlim([0 1000]);
grid on;

settling_time_value = find(abs(y_step - y_step(end)) >= 0.01, 1, 'last') + 1;
fprintf("Setteling sample (H2): %d\n", settling_time_value);
%% 3.2.c
f = [4, 8 ,12];
fs = 400;

t1 = 0:(1/fs):(2-1/fs);
t2 = 2:(1/fs):(4-1/fs);
t3 = 4:(1/fs):(6-1/fs);
x = [cos(2*pi*f(1)*t1),cos(2*pi*f(2)*t2),cos(2*pi*f(3)*t3)];

% H1:
figure('Name', 'H1');
plot(x, "LineWidth", 1.5);
title("Filtered Output (H1)");
xlabel("n");
ylabel("Amplitude");
grid on;
hold on;

y1 = filter(b1, a1, x);
plot(y1, "LineWidth", 1.5);
hold off;
legend('x[tn]', 'y[tn]')

% H2:
figure('Name', 'H2');
plot(x, "LineWidth", 1.5);
title("Filtered Output (H2)");
xlabel("n");
ylabel("Amplitude");
grid on;
hold on;

y2 = filter(b2, a2, x);
plot(y2, "LineWidth", 1.5);
hold off;
legend('x[tn]', 'y[tn]')

%% 3.2.d
N = 300;
n = 0:N;

impr1 = impz(b1,a1,n);
impr2 = impz(b2,a2,n);

[max_h1, indx] = max(impr1);

[max_h2, indx] = max(impr2);
fprintf("Max H1: %f\nMax H2: %f\n", max_h1, max_h2)

figure('name', "Impulse Res H1")
plot(n, abs(H1), 'LineWidth', 1.5);
grid on;
xlabel('time');
ylabel('|H1(z)|');
title('Impulse Response H1');


figure('name', "Impulse Res H2")
plot(n, abs(H2), 'LineWidth', 1.5);
grid on;
xlabel('time');
ylabel('|H2(z)|');
title('Impulse Response H2');

%% 3.2.f
fs  = 400;
[h1, w1] = freqz(b1,a1,2048);
[h2, w2] = freqz(b2,a2,2048);

% Plotting frequency response of H1 and H2
figure('name', "Notch Filter Responses")
plot(w1/pi*fs/2, abs(h1), 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H(f)|'); 
hold on;
title('Frequency Responses of H1 and H2');
plot(w2/pi*fs/2, abs(h2), '--', 'LineWidth', 1.5);
legend('H1(f) Δf = 4', 'H2(f) Δf = 0.5', 'Location','best');