clc
clear
close all

f_c = 4.3e9;
c = 3e8;
lambda=c/f_c;
f_lowIF = 5e4;                  %Low intermediate frequency

fs = 5e5;                       %sampling frequency
T = 10;
Ns = fs*T;
n = 0:Ns-1;
t = 1/fs .* n;
f_axis = n*1/T;
Delta_f = 1/T;

A_heart = 1e-3;                 %Amplitude of heart beats
f_heart = 1.25;
A_resp =  0e-3;             %Amplitude of respiration
f_resp = 0;

h_t = A_heart*sin(2*pi*f_heart*t);  %signal of heart beat
r_t = A_resp*sin(2*pi*f_resp*t);    %signal of respiration

IF_t = exp(1i*(2*pi*f_lowIF.*t-(h_t+r_t)*4*pi/lambda+ 0.3*wgn(1,Ns,1))) ;

N_segment = 1000;
N_IF = N_segment/10+1;
N_piece = Ns/N_segment;

segment_matrix = reshape(IF_t,[N_segment,Ns/N_segment]);
sm_fft = fft(segment_matrix);

phi = angle(sm_fft(N_IF,:));

for n=1:length(phi)-1
    if phi(n+1) - phi(n) > pi/2
        phi(n+1) = phi(n+1) - pi;
    elseif phi(n+1) - phi(n) < -pi/2
        phi(n+1) = phi(n+1) + pi;
    end
end

phi = detrend(phi,0);
trajectory = phi * lambda/4/pi;

D_t = T/N_piece;

D_f = 1/D_t;
f_axis = (0:N_segment-1)/D_t;

t_axis = (0:N_piece-1)/D_f;

Spectrum = abs(fft(phi)); 

%%
figure(1)
mesh(t_axis,f_axis,abs(sm_fft));
colorbar;
view(0,90);
title('Low-IF Signal of CW Radar')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%%
figure(2)
ptitle = ['N = ',num2str(N_segment)];
plot(t_axis,trajectory);
xlabel('Time (s)')
ylabel('motion (m)')
title(ptitle)



