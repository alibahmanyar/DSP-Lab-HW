%%
close all;
clear;
clc;
%% part 1
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

%% part 2
% part 1.2
r = rand(1, length(y)) - 0.5; % Noise
noised_signal = y + r;        % Add noise to sine, y defined in section 1.1

% Plotting signals
figure('Name', 'Normal Sine signal VS Noisy Sine Signal', 'Position', [300 200 700 500]);
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

%% part 3
% part 1.3:
% A moving average filter introduces a delay because it computes the average of past data points, which inherently involves looking back in time, causing the output to lag behind the input signal.

M1 = 0;
M2 = 20;
len = M2 + M1 + 1; % Window length

window = ones(1, len) / len; % Moving average window

y2 = conv(noised_signal, window);   % Noised signal defined in section 1.2, Line 30
t2 = 0:0.01:(0.01 * (length(y2) - 1));

% Plotting section 1.3.1
figure('Name', 'Moving Avg signals', 'Position', [300 200 700 500]);
hold on;
plot(t,noised_signal);              % Noise signal
plot(t2, y2);                       % Moving Avg M1=0, M2=20

% part 1.3.2:
M1 = 10;
M2 = 10;
len = M2 + M1 + 1; % Window length

window = ones(1, len) / len; % Moving average window

y3 = conv(noised_signal(M1:end), window);
t3 = 0:0.01:(0.01 * (length(y3) - 1));

plot(t3, y3);                       % Moving Avg M1=0, M2=20
hold off;
title('Three signals together'); xlabel('Time (Second)'); ylabel('Signal Amplitude');
legend("Noisy Signal", "Moving Avg M1=0, M2=20", "Moving Avg M1=10, M2=10");

%% part 4
% part 1.4.1
M1 = 0;
M2 = 20;
len = M2 + M1 + 1; % Window length

b = ones(1,len)/len;
a = 1;
y4_1 = filter(b, a, noised_signal);
t4_1 = 0:0.01:(0.01 * (length(y4_1) - 1));

% Plotting section 1.4
figure('Name', 'Moving Avg signals by using filter function', 'Position', [300 200 700 500]);
hold on;
plot(t,noised_signal);              % Noise signal
plot(t4_1, y4_1);                   % filtered signal
plot(t2,y2);                        % Plotting convolved signal M1=0, M2=20
plot(t3,y3);                        % Plotting convolved signal M1=10, M2=10
title('All signals together'); xlabel('Time (Second)'); ylabel('Signal Amplitude');
legend("Noisy Signal", "filtered signal", "Moving Avg M1=0, M2=20", "Moving Avg M1=10, M2=10");
hold off;

% Plotting section 1.3
figure('Name', 'plotting in seperated subplots', 'Position', [300 200 700 500]);
subplot(2,2,1);               % Plotting noisy signal
plot(t,noised_signal)
title('Noisy signal');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(2,2,2);               % Plotting convolved signal M1=0, M2=20
plot(t2,y2)
title('convolved signal M1=0, M2=20');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(2,2,3);               % Plotting convolved signal M1=10, M2=10
plot(t3,y3)
title('convolved signal M1=10, M2=10');xlabel('Time (Second)');ylabel('Signal Amplitude');
subplot(2,2,4);               % filtered signal
plot(t4_1,y4_1)
title('output signal of filter function');xlabel('Time (Second)');ylabel('Signal Amplitude');

%% part 5
% part 1.5
w0 = pi/32; n = 100;
[y5,t5] = singen(w0,n);     %singen function is defined in singen.m
figure('Name', 'singen');
stem(t5,y5, '.')
title('singen function output');xlabel('Time (Second)');ylabel('Signal Amplitude');


%% part 6
close all;
clear;
clc;

figure('Name', 'cos(2*pi*t6) + cos(8*pi*t6) + cos(12*pi*t6)');
t6 = 0:0.01:4;
y6 = cos(2*pi*t6) + cos(8*pi*t6) + cos(12*pi*t6); %Original Signal
plot(t6,y6);
title('cos(2*pi*t6) + cos(8*pi*t6) + cos(12*pi*t6)');xlabel('Time (Milisecond)');ylabel('Signal Amplitude');
hold on;

t62 = 0:0.2:4; % 5khz => 0.0002s = 0.2ms
y62 = cos(2*pi*t62) + cos(8*pi*t62) + cos(12*pi*t62); %Same signal, 5Khz sampled time vector
stem(t62,y62);

%figure()
%y8 = idealfilter()
%t8 = 0:0.01:(0.01*(length(y8)-1));

%plot(t8,y8)

%% part 7
close all;
clear;
clc;

fs = 100;
t7 = -5:1/fs:5;
y7 = sinc(5*t7) .^ 2;

% Spectrum of original singal
FT_x_0 = (1 / fs) * abs(fftshift(fft(y7)));
f_axis_0 = linspace(-fs / 2, fs / 2, length(FT_x_0));


for sr=[4 5 10 20] % sr: sampling rate
    figure('Name', 'Spectrum of original vs sampled signal');
    
    % Keep every fs/sr sample, other samples are set to 0
    % This way, the length of y will remain the same
    y7_2 = zeros(1, length(y7));
    y7_2(1:fs/sr:end) = y7(1:fs/sr:end); 


    % Spectrum of sampled singal
    FT_x_1 = (1 / fs) * abs(fftshift(fft(y7_2)));
    f_axis_1 = linspace(-fs / 2, fs / 2, length(FT_x_1));

    subplot(2,1,1)
    plot(f_axis_0, abs(FT_x_0), 'LineWidth', 1.5);
    title('Spectrum of original signal');xlabel('Frequency (Hz)');ylabel('Amplitude');

    subplot(2,1,2)
    plot(f_axis_1, abs(FT_x_1), 'LineWidth', 1.5);
    title(sprintf('Spectrum of sampled signal sr=%d', sr));xlabel('Frequency (Hz)');ylabel('Amplitude');
end


%% part 9
close all;
clear;
clc;

% Constructing the signal
f = [pi/16 5*pi/16 9*pi/16 13*pi/16];

fs = 1e3;
t9 = 1:1/fs:100;
y9 = zeros(1, length(t9));

for f1=f
    y9 = y9 + cos(2*f1*t9);
end


% Plottign the spectrum of the original signal
figure('Name', 'Spectrum of original signal');
FT_y9 = (1/fs) * abs(fftshift(fft(y9)));
f_axis = linspace(-fs / 2, fs / 2, length(FT_y9));
plot(f_axis,FT_y9, 'LineWidth', 2)
xlim([-1,1])
title('Spectrum'); xlabel('Freq');ylabel('Amplitude');


% Loading filters coefficients from xls file
filter_coef_1 = readmatrix('filters.xls', 'Sheet', 1);
filter_coef_2 = readmatrix('filters.xls', 'Sheet', 2);



% Plotting amplitude of frequency responses
figure('Name', 'Frequency response of filters');

for i = 1:4
    subplot(4,2, 2*i-1)
    
    [h, w] = freqz(filter_coef_1(i, :));
    w = w ./ pi;

    plot(w, abs(h));
    title(sprintf('Frequency response of analysis filter %d', i));
    xlabel('Normalized Frequency');ylabel('Amplitude');
end

for i = 1:4
    subplot(4,2,2*i)
    
    [h, w] = freqz(filter_coef_2(i, :));
    w = w ./ pi;
    
    plot(w, abs(h));
    title(sprintf('Frequency response of synthesis filter %d', i));
    xlabel('Normalized Frequency');ylabel('Amplitude');
end

% Plotting phase of frequency responses
figure('Name', 'Phase response of filters');

for i = 1:4
    subplot(4,2, 2*i-1)
    
    [h, w] = freqz(filter_coef_1(i, :));
    w = w ./ pi;
    
    plot(w, (angle(h)));
    title(sprintf('Frequency response of analysis filter %d', i));
    xlabel('Normalized Frequency');ylabel('Phase');
end

for i = 1:4
    subplot(4,2,2*i)
    
    [h, w] = freqz(filter_coef_2(i, :));
    w = w ./ pi;
    
    plot(w, (angle(h)));
    title(sprintf('Frequency response of synthesis filter %d', i));
    xlabel('Normalized Frequency');ylabel('Phase');
end