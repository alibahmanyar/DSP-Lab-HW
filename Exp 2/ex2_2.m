%%
close all;
clear;
clc;

%% 2.3.a
% Loading the audio file
[x, Fs] = audioread('Audio01.wav');
T = 10;
t = 0: 1/Fs:T-1/Fs;

% Constructing the FFT spectrum of the signal
L = length(x);
f_x = (Fs/L) * (-L/2:L/2-1);
fft_x = fftshift(fft(x))/L;

% Plotting signal in time and freq domain
figure('Name', 'Original audio signal');
subplot(2, 1, 1);
plot(t, x, 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Signal in Time Domain');
subplot(2,1,2) ;
plot(f_x, abs(fft_x), 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('Amplitude');
title('Signal in Frequency Domain');

% Playing the audio
sound(x, Fs);
pause(length(x) * (1/ Fs));

%% 2.3.b
load('filter.mat');

freqz(Num); % Visualising filter freq-phase response
y_0 = filter(Num,1,x);

% Constructing the FFT spectrum of the signal
L_y_0 = length(y_0);
f_y_0 = (Fs/L_y_0) * (-L_y_0/2:L_y_0/2-1);
fft_y_0 = fftshift(fft(y_0))/L_y_0;

% Plotting signal in time and freq domain
figure('Name', 'Filtered audio signal');
subplot(2, 1, 1);
plot(t, y_0, 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Signal in Time Domain');
subplot(2,1,2) ;
plot(f_y_0, abs(fft_y_0), 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('Amplitude');
title('Signal in Frequency Domain');


%% 2.3.c
f0 = 10000;
n   = 1:length(x);
s_n = 2*cos(2*pi*f0*n/Fs);

y_1 = y_0 .* s_n';
y_2 = filter(Num,1,y_1);

% Constructing the FFT spectrum of the signal
L_y_2 = length(y_2);
f_y_2 = (Fs/L_y_2) * (-L_y_2/2:L_y_2/2-1);
fft_y_2 = fftshift(fft(y_2))/L_y_2;

% Plotting signal in time and freq domain
figure('Name', 'Modulated audio signal');
subplot(2, 1, 1);
plot(t, y_2, 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Signal in Time Domain');
subplot(2,1,2) ;
plot(f_y_2, abs(fft_y_2), 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('Amplitude');
title('Signal in Frequency Domain');


% Playing the audio
sound(y_2, Fs);
pause(length(y_2) * (1/ Fs));


%% 2.3.d
y_3 = filter(Num,1,y_2);
y_4 = y_3 .* s_n';
y_5 = filter(Num,1,y_4);

% Constructing the FFT spectrum of the signal
L_y_5 = length(y_5);
f_y_5= (Fs/L_y_5) * (-L_y_5/2:L_y_5/2-1);
fft_y_5 = fftshift(fft(y_5))/L_y_5;

% Plotting signal in time and freq domain
figure('Name', 'Reconstructed audio signal');
subplot(2, 1, 1);
plot(t, y_5, 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Signal in Time Domain');
subplot(2,1,2) ;
plot(f_y_2, abs(fft_y_5), 'LineWidth',1.5) ;
grid on;
xlabel('freq (Hz)');
ylabel('Amplitude');
title('Signal in Frequency Domain');

% Playing the audio
sound(y_5, Fs);

% Calculating MSE and MAE
MAE = mean(abs(x - y_5))
MSE = mean((x - y_5) .^ 2)
