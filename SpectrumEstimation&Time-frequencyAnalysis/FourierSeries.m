clc
clear all
close all

Fs = 1e3;
Ns = 1e3;

PA = Param(Fs,Ns);

f = 1;
x = square(2*pi*f*PA.t);

c = [];
N = 200;        % Change this para to set the number of harmonic

for i=1:N
    c(i) = 2/Fs * sum(x .* exp(-1i*i*2*pi*f*PA.t));
end

x_t = 0;
for i=1:N
    x_t =  x_t + c(i)*exp(1i*i*2*pi*f*PA.t);
end

ptitle = ['Gibbs Phenomenon N = ',num2str(N)];

figure;
plot(PA.t, x)
title(ptitle)
xlabel('Time(s)')
ylabel('Amplitude')
hold on
plot(PA.t, real(x_t))
hold off
legend('Origin Signal',' Composition of harmonic')

