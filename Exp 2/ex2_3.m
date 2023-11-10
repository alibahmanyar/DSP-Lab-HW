%% 2.5
loc = linspace(0,1,2^10);
[x,x_noisy] = wnoise('doppler',10,7);

figure('Name', 'Clean Vs. Noisy Doppler');

subplot(2,1,1);
plot(loc,x, 'g', 'LineWidth',1.5); 
grid on; title('Clean Doppler');

subplot(2,1,2);
plot(loc,x_noisy, 'c', 'LineWidth',1.5);
grid on; title('Noisy Doppler');

[cA1,cD1] = dwt(x_noisy,'db1');
loc = linspace(0,1,length(cA1));
