clc
clear all
close all

Fs = 1e2;
Ns = 1.5e2;
PA = Param(Fs,Ns);

f = 1;
x0 = 0*PA.t-1;
x1 = cos(2*pi*f*PA.t);



s0 = [x0 x1 x0];
s1 = [x1 x1 x1];

figure;
plot(s0)
title('Signal In Time Window')
ylim([-1 1.1])

figure;
plot(s1)
title('Signal After Period Extension')
ylim([-1 1.1])

%%
figure;
plot(PA.f,abs(fft(x1))/Ns);
title('Spectrum')
xlabel('Frequency (Hz)')
ylabel("Amplitude")
ylim([0 1.1])

%% rectangular windows

rectwin = [zeros(1,Ns) ones(1,Ns) zeros(1,Ns)];
figure;
plot(rectwin)
ylim([0 1.1])
title('Rectangular Windows')

%% Fence effect

Fs = 1e2;
Ns = 1e2;
PA = Param(Fs,Ns);

delta_f = 1/(Ns/Fs);
f0 = 9.5;
x_t = exp(1i*2*pi*f0*PA.t);

figure;
stem(PA.f,abs(fft(x_t))/Ns)
ptitle = ['Delta frequency = ',num2str(delta_f), ', signal frequency =',num2str(f0)];
title(ptitle)
ylim([0 1.1])

