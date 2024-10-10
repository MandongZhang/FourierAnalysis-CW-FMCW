close all
clc
clear all
% Signal Parameter
Fs = 2e3;
T = 0.1;
Ns = Fs*T;
PA = Param(Fs,Ns);
t = PA.t;
f = PA.f;
n = PA.n;

Delta_f = 1/T;
fc = 1e2;

% Complex Signal
x_t = exp(1i*(2*pi*fc*t+pi/4));

X_f = fft(x_t)/Ns;
X_A = abs(X_f);
X_Phi = atan2(imag(X_f),real(X_f));

figure;
plot(f,X_A)
title("Amplitude of Signal")
Index_max = find(X_A==max(X_A));
text( f(Index_max), X_A(Index_max),['(f=',num2str(f(Index_max)),')']);
ylim([0 1.1])

figure;
plot(f,X_Phi)
title('Origin Phase Spectrum')
Index_max = find(X_Phi==max(X_Phi));


%% Phase Spectrum Adjustment

% Due to the bit precision of double type in calculation, some frequency
% points exist small which will be considered in the phase calculation.
% The Scheme is to set a threshold of Amplitude Spectrum, and the amplitude
% lower than threshold is set to 0

X_f2 = X_f;
threshold = max(X_A)/10000;
X_f2(X_A<threshold) = 0;
X_Phi2 = atan2(imag(X_f2),real(X_f2))/pi;
figure;
plot(f,X_Phi2);
title('Adjusted Phase Spectrum');
xlabel('Frequency(Hz)')
ylabel('Ampltide(Pi)')



%% Charactistic of Frequency Shift in FFT
f_shift = 101;
x_tshift = x_t.*exp(1i*2*pi*f_shift.*t);    %frequency shift part should  multiply with time-domain signal 
X_fshift = fft(x_tshift);
X_Ashift = abs(X_fshift);

% Phase Spectrum filter
X_fshift2 = X_fshift;
threshold = max(X_Ashift)/1000;
X_fshift2(abs(X_fshift)<threshold)=0;

X_Phishift = atan2(imag(X_fshift2), real(X_fshift2))/pi;

figure
plot(f,abs(X_fshift));
xlabel('Frequency(Hz)')
ylabel('Amplitude');

figure;
plot(f,X_Phishift)
title('Phase Spectrum')

%%






