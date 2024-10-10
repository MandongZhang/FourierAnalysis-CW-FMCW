clc
clear all
close all

Fs = 2e2;
Ns = 5e2;
PA = Param(Fs,Ns);

f1 = 10;
f2 = 20;
f3 = 50;

x1 = exp(1i*2*pi*f1*PA.t);
x2 = exp(1i*2*pi*f2*PA.t);
x3 = exp(1i*2*pi*f3*PA.t);


x_sum1 = [x1 x2 x3];
x_sum2 = [x2 x3 x1];

Ns_sum = 3*Ns;
PA_sum = Param(Fs,Ns_sum);

figure;
plot(PA_sum.t, x_sum1);
xlabel('Time (s)')
title('Signal 1')

figure;
plot(PA_sum.t, x_sum2);
xlabel('Time (s)')
title('Signal 2')

%%
figure;
plot(PA_sum.f, abs(fft(x_sum2))/Ns);
title('Spectrum of all data')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% STFT 
% window = hamming(256);
% noverlap = 128;

[stft1, f, t_stft] = stft(x_sum1,Fs,'FFTLength',128);

figure;
surf(t_stft, f, abs(stft1), 'EdgeColor', 'none');
axis tight;
ylim([-100 100])
% view(0, 90);
xlabel('时间 (s)');
ylabel('频率 (Hz)');
title('STFT 幅度谱');
colorbar;



