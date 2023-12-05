close all;
clear;
clc;

%%
% H3
b1 = [0.030469 0 -0.030469];
a1 = [1 -1.923772 0.939063];
% H4:
b2 = [0.003912 0 -0.003912];
a2 = [1 -1.976468 0.992177];

fs = 400;
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
subplot(2,1,1);
plot(x, "LineWidth", 1.5);
title("Input signal");
xlabel("n");
ylabel("Amplitude");
grid on;

y1 = filter(b1, a1, x);
subplot(2,1,2);
plot(y1, "LineWidth", 1.5);
title("Filtered Output (H1)");
xlabel("n");
ylabel("Amplitude");
grid on;

% H2:
figure('Name', 'H2');
subplot(2,1,1);
plot(x, "LineWidth", 1.5);
title("Input signal");
xlabel("n");
ylabel("Amplitude");
grid on;

y2 = filter(b2, a2, x);
subplot(2,1,2);
plot(y2, "LineWidth", 1.5);
title("Filtered Output (H2)");
xlabel("n");
ylabel("Amplitude");
grid on;

%% 3.2.d
N = 300;
n = 0:N;

impr1 = impz(b1,a1,n);
impr2 = impz(b2,a2,n);

[max_h1, indx] = max(impr1);

[max_h2, indx] = max(impr2);
fprintf("Max H1: %f\nMax H2: %f\n", max_h1, max_h2)

figure('name', "Impulse Res H1")
plot(n, abs(impr1), 'LineWidth', 1.5);
grid on;
xlabel('time');
ylabel('|H1(z)|');
title('Impulse Response H1');


figure('name', "Impulse Res H2")
plot(n, abs(impr2), 'LineWidth', 1.5);
grid on;
xlabel('time');
ylabel('|H2(z)|');
title('Impulse Response H2');

%% 3.2.f-g
fs  = 400;
Nf = 4096;
[h1, w1] = freqz(b1, a1, Nf);
[h2, w2] = freqz(b2, a2, Nf);

% Plotting frequency response of H1 and H2
figure('name', "Peaking Filter Responses")
plot(w1/pi*fs/2, abs(h1), 'LineWidth', 1.5);
xlim([0, 20]);
grid on;
xlabel('freq'); 
ylabel('|H(f)|'); 
hold on;
title('Frequency Responses of H1 and H2');
plot(w2/pi*fs/2, abs(h2), '--', 'LineWidth', 1.5);

% Finding section of array which resides in the bandwidth of response
bw1 = find(abs(abs(h1) - .5^.5*max(abs(h1))) < 0.02);
bw2 = find(abs(abs(h2) - .5^.5*max(abs(h2))) < 0.02);

fli1 = bw1(1); % F_l_1 index
fhi1 = bw1(end); % F_h_1 index

fli2 = bw2(1); % F_l_2 index
fhi2 = bw2(end); % F_h_2 index

plot(w1(fli1)/pi*fs/2, abs(h1(fli1)), 'ro', 'LineWidth', 1.5);
plot(w1(fhi1)/pi*fs/2, abs(h1(fhi1)), 'ro', 'LineWidth', 1.5);

plot(w2(fli2)/pi*fs/2, abs(h2(fli2)), 'ko', 'LineWidth', 1.5);
plot(w2(fhi2)/pi*fs/2, abs(h2(fhi2)), 'ko', 'LineWidth', 1.5);


legend('H1(f) Δf = 4', 'H2(f) Δf = 0.5','fL_1','fH_1','fL_2','fH_2', 'Location','best');

fprintf("H1:\nFh-Fl=%f\n\nH2:\nFh-Fl=%f\n", (fhi1-fli1)*(fs/Nf/2), (fhi2-fli2)*(fs/Nf/2))