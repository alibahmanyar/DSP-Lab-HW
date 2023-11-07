%%
close all;
clear;
clc;

%% 2.3.a
% Loading the audio file
[x, Fs] = audioread('Audio01.wav');

% Plotting signal in time and freq domain
plot_time_freq(x, Fs, 'Original audio signal', 'Signal in Time Domain', 'Signal in Frequency Domain');

% Playing the audio
sound(x, Fs);
pause(length(x) * (1/ Fs));

%% 2.3.b
load('filter.mat');

% Plotting filter in time domain
figure('Name', 'Time Domain Representation of Filter');
stem(Num, 'LineWidth', 1.5);
grid on;
xlabel('Time(Sec)') ;
ylabel('Amplitude') ;
title('Time Domain Representation of Filter');

figure('Name', 'Filter freq-phase response');
freqz(Num); % Visualising filter freq-phase response

y_0 = filter(Num,1,x); % Applying filter


% Plotting signal in time and freq domain
plot_time_freq(y_0, Fs, 'Filtered audio signal', 'Filtered Signal in Time Domain', 'Filtered Signal in Frequency Domain');
%% 2.3.c
f0 = 10000;
n   = 1:length(x);
s_n = 2*cos(2*pi*f0*n/Fs);

y_1 = y_0 .* s_n';
y_2 = filter(Num,1,y_1);

% Plotting signals in time and freq domain
plot_time_freq(y_1, Fs, 'Signal Multiplied by Carrier', 'Signal Multiplied by Carrier in Time Domain', 'Signal Multiplied by Carrier in Frequency Domain');
plot_time_freq(y_2, Fs, 'Modulated Audio Signal', 'Modulated Audio Signal in Time Domain', 'Modulated Audio Signal in Frequency Domain');

% Playing the audio
sound(y_2, Fs);
pause(length(y_2) * (1/ Fs));

%% 2.3.d
y_3 = filter(Num,1,y_2);
y_4 = y_3 .* s_n';
y_5 = filter(Num,1,y_4);

% Plotting signals in time and freq domain
plot_time_freq(y_4, Fs, 'Scrmabled Signal Multiplied by Carrier', 'Scrmabled Signal Multiplied by Carrier in Time Domain', 'Scrmabled Signal Multiplied by Carrier in Frequency Domain');
plot_time_freq(y_5, Fs, 'Descrmabled Audio Signal', 'Descrmabled Audio Signal in Time Domain', 'Descrmabled Audio Signal in Frequency Domain');

% Playing the audio
sound(y_5, Fs);

% Calculating MSE and MAE
MAE = mean(abs(x - y_5))
MSE = mean((x - y_5) .^ 2)


function plot_time_freq(y, Fs, title1, title2, title3)
    T = length(y) * (1/ Fs);
    t = 0: 1/Fs:T-1/Fs;
    

    % Constructing the FFT spectrum of the signal
    L_y = length(y);
    f_y = (Fs/L_y) * (-L_y/2:L_y/2-1);
    fft_y = fftshift(fft(y))/L_y;
    
    % Plotting signal in time and freq domain
    figure('Name', title1);
    subplot(2, 1, 1);
    plot(t, y, 'LineWidth', 1.5);
    grid on;
    xlabel('Time(Sec)') ;
    ylabel('Amplitude') ;
    title(title2);
    subplot(2,1,2) ;
    plot(f_y, abs(fft_y), 'LineWidth',1.5) ;
    grid on;
    xlabel('freq (Hz)');
    ylabel('Amplitude');
    title(title3);
end