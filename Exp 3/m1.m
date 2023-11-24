clear all;
close all;
clc;

%% 3-1-a:
clear all;
close all;
clc;

n = 500;
f0 = 500;
fs = 10^4;

w = linspace(0, pi, n);
R = [0.8, 0.9, 0.99];

w0 = f0 * 2 * pi / fs;

for r=R
    G = (1-r)*(1-2*r*cos(w0)+r^2)^.5;
    b = [G];
    a = [1, -2*r*cos(2*w0), r^2];
    
    figure('Name', sprintf("R=%d", r))
    freqz(b, a)
    % [Amp, physconst] = freqz(b, a, n)
    % plot(w, Amp);
end


%% 3-1-b:
clear all;
close all;
clc;

n = 300;
f0 = 500;
fs = 10^4;

w = linspace(0, pi, n);
R = [0.8, 0.9, 0.99];

w0 = f0 * 2 * pi / fs;
y = zeros(1,n);

for r=R

    G = (1-r)*(1-2*r*cos(w0)+r^2)^.5;
    b = [G];
    a = [1, -2*r*cos(2*w0), r^2];
    
    y(1) = G;
    y(2) = -a(2)*y(1);
    

    hn = zeros(1, n);
    for k = 1:n
        hn(k) = (G/(sin(w0)))*(r^k)*sin(w0*(k + 1));
    end
    

    for k = 3:n
        y(k) = -a(2)*y(k-1) - a(3)*y(k-2);
    end

    
    
    figure('Name', sprintf("R=%d", r))
    plot(y);
    hold on;
    plot(hn);
    hold off;
end

%% 3-1-c:
clear all;
close all;
clc;

f0 = 500; % Hz
fs = 10000; % Hz
w  = 2 * pi * f0 / fs;
R = [0.80, 0.90, 0.99];

N = 300;
n  = 0:1:N;
s  = cos (w*n);
v  = randn(1,length(n));
x  = s + v;

figure('Name', 'Noisy Signal');
plot(n, x, 'r');
grid on;
xlabel('n');
ylabel('x[n]');
title('Noisy Signal');

for ll = 1:length(R)
    G = (1 - R(ll)) * (1 - 2*R(ll)*cos(2*w) + R(ll)^2)^0.5;
    b = G;
    a = [1, -2*R(ll)*cos(w), R(ll)^2];
    y = zeros(1,N);
    for mm = 1:N
        if (mm == 1)
            y(mm) = G*x(mm) ;
            w1 = y(1) ;
        elseif(mm == 2)
            y(mm) = -a(2)*w1 + G*x(mm) ;
            w2 = w1;
            w1 = y(2) ;
        else
            y(mm) = -a(2)*w1 - a(3)*w2 + G*x(mm) ;
            w2 = w1;
            w1 = y(mm);
        end
    end
    figure(8)
    subplot(length(R), 1, ll);
    plot(n, s, 'g');
    hold on;
    plot(n(1:end-1), y, 'r');
    grid on;
    xlabel('n');
    ylabel('y[n]');
    title("R = " + num2str(R(ll)));
end
%% 3-1-d:
clear all;
close all;
clc;

f0 = 500; % Hz
fs = 10000; % Hz
w  = 2 * pi * f0 / fs;
R = [0.80, 0.90, 0.99];
N = 300;
n  = 0:1:N;
s  = cos (w*n);
v  = randn(1,length(n));
x  = s + v;

for ll = 1:length(R)
    G = (1 - R(ll)) * (1 - 2*R(ll)*cos(2*w) + R(ll)^2)^0.5;
    b = G;
    a = [1, -2*R(ll)*cos(w), R(ll)^2];
    y_v = zeros(1,N);
    for mm = 1:N
        if (mm == 1)
            y_v(mm) = G*v(mm) ;
            w1 = y_v(1) ;
        elseif(mm == 2)
            y_v(mm) = -a(2)*w1 + G*v(mm) ;
            w2 = w1;
            w1 = y_v(2) ;
        else
            y_v(mm) = -a(2)*w1 - a(3)*w2 + G*v(mm) ;
            w2 = w1;
            w1 = y_v(mm);
        end
    end
    if(ll == 1)
        figure(9)
        subplot(length(R)+1, 1, ll);
        plot(n, v, 'b');
        hold on;
        grid on;
        xlabel('n');
        ylabel('v[n]');
        title('noise');
        legend('v[n]');
    end
    figure(9)
    subplot(length(R)+1, 1, ll+1);
    plot(n(1:end-1), y_v, 'c');
    grid on;
    xlabel('n');
    ylabel('y_v[n]');
    title("filtered noise, R = " + num2str(R(ll)));
    legend('y_v[n]');
end
%% 3-1-e:

clear all;
close all;
clc;

f0 = 500; % Hz
fs = 10000; % Hz
w  = 2 * pi * f0 / fs;
R = [0.80, 0.90, 0.99];
N = 300;
n  = 0:1:N;
s  = cos (w*n);
v  = randn(1,length(n));
x  = s + v;
for ll = 1:length(R)
    G = (1 - R(ll)) * (1 - 2*R(ll)*cos(2*w) + R(ll)^2)^0.5;
    b = G;
    a = [1, -2*R(ll)*cos(w), R(ll)^2];
    y_v = zeros(1,N);
    for mm = 1:N
        if (mm == 1)
            y_v(mm) = G*v(mm) ;
            w1 = y_v(1) ;
        elseif(mm == 2)
            y_v(mm) = -a(2)*w1 + G*v(mm) ;
            w2 = w1;
            w1 = y_v(2) ;
        else
            y_v(mm) = -a(2)*w1 - a(3)*w2 + G*v(mm) ;
            w2 = w1;
            w1 = y_v(mm);
        end
    end
    NRR = std(y_v)/std(v);
    NRR_th = (1+R(ll)^2)/((1+R(ll))*(1+2*R(ll)*cos(w)+R(ll)^2));

    figure(10);
    stem(R(ll), NRR, 'r');
    hold on;
    stem(R(ll), NRR_th, 'g');
    grid on;
    xlabel('R');
    ylabel('NRR');
    title('NRR from eq. Vs. Theory');
    legend('NRR', 'NRR Theory');
end
