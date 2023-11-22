clear all;
close all;
clc;

%%
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


%%
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
