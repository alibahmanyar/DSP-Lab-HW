%% 2.5
loc = linspace(0,1,2^10);
[x,x_noisy] = wnoise('doppler',10,7);

figure('Name', 'Clean Vs. Noisy Doppler');

subplot(2,1,1);
plot(loc,x, 'g', 'LineWidth',1.5); 
grid on; title('Clean Doppler');
xlabel('Time(Sec)') ;
ylabel('Amplitude');

subplot(2,1,2);
plot(loc,x_noisy, 'c', 'LineWidth',1.5);
grid on; title('Noisy Doppler');
xlabel('Time(Sec)') ;
ylabel('Amplitude');

[cA1,cD1] = dwt(x_noisy,'db1');
loc = linspace(0,1,length(cA1));

figure('Name', '3 Levels');
subplot(3,1,1);
plot(loc,cA1, 'r', 'LineWidth',1.5);
grid on; title('Level 1');
xlabel('Time(Sec)') ;
ylabel('Amplitude');

[cA2,cD2] = dwt(cA1,'db1');
loc = linspace(0,1,length(cA2));
subplot(3,1,2);
plot(loc,cA2, 'm', 'LineWidth',1.5);
grid on; title('Level 2');
xlabel('Time(Sec)') ;
ylabel('Amplitude');

[cA3,cD3] = dwt(cA2,'db1');
loc = linspace(0,1,length(cA3));
subplot(3,1,3);
plot(loc,cA3, 'g', 'LineWidth',1.5);
grid on; title('Level 3');
xlabel('Time(Sec)') ;
ylabel('Amplitude');