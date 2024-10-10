clc
clear all
close all

rng(42);

c = 3e8;
fc = 5.8e9;
BW = 1.6e8;
Fs = 8e3;
Ns_chirp = 40;
Nc = 2400;
Fchirp = 200;
Tchirp = 1/Fchirp;
T =  Nc*Tchirp;
Ns = T*Fs;
PA = Param(Fs,Ns);
t = PA.t;
f = PA.f;
t_fast = (0:Ns_chirp-1)/Fs;
f_fast =(0:Ns_chirp-1)/Ns_chirp*Fs;

lambda = c/fc;
delta_range = c/2/BW;
axis_range = (0:Ns_chirp-1)*delta_range;
delta_f = Fs/Ns_chirp;
delta_phi_FT = 4*pi*fc*delta_range/c;

%% Moving target modeling
R  = 5*delta_range-5/10*delta_range;
fm = 0.5;
Am = 1e-2;
x_t = Am*(cos(2*pi*fm*t));
R_t = R + x_t; 
tau = 2*R_t/c;

v = 1;
C_t = 1+v*t;
tau_v = 2*C_t/c;

%% Tx Rx IF signal generation


f_slow = (0:Nc-1)/Tchirp;
t_slow = (0:Nc-1)*Tchirp;
Slope = BW / Tchirp;

Tx = exp(1i*(2*pi*fc.*(t)+pi*Slope.*mod(t, Tchirp).^2));
Rx1 = exp(1i*(2*pi*fc.*(t-tau) + Slope*pi.*mod((t-tau), Tchirp).^2) );
Rx2 =  exp(1i*(2*pi*fc.*(t-tau_v) + Slope*pi.*mod((t-tau_v), Tchirp).^2) );

% Sb_t = 1*Tx ./ Rx1+ 0*Tx./Rx2 + wgn(1,Ns,3);  
Sb_t = 1*Tx ./ Rx1+ 0*Tx./Rx2 + wgn(1,Ns,2) ;
% Sb_t = awgn(Sb_t,-10);

Sb_f = abs(fft(Sb_t));
%% reshape matrix and add win
ChirpArray = reshape(Sb_t,Ns_chirp,Nc)';
% Win = blackman(Ns_chirp)';
% ChirpArray = ChirpArray.*Win;

%% Range FFT
FT_fft =fft(ChirpArray')';
figure;
mesh(axis_range,t_slow,abs(FT_fft))
axis tight;
view(0,90);
colorbar;
xlabel('Range')
ylabel('Slow Time Index')
title('Range-FFT')

%% Amplitude Spectrum and Phase Spectrum of A signal chirp
ChirpSlice = 1;
Phi_SignalChirp = unwrap(atan2(imag(FT_fft(ChirpSlice,:)),real(FT_fft(ChirpSlice,:))));
% threshold = max(abs(FT_fft(ChirpSlice,:)))/1000;
% Phi_SignalChirp(abs(FT_fft(ChirpSlice,:))<threshold ) = 0;

figure;
subplot(211)
plot(axis_range,abs(FT_fft(ChirpSlice,:)))
xlabel('Range')
title('Amplitude Spectrum')
subplot(212)
plot(axis_range,Phi_SignalChirp)
xlabel('Range')
title('Phase Spectrum')

%% CI (Coherent Integration)
FT_fft_CI = sum((FT_fft))/Ns_chirp/Nc;
figure;
plot(axis_range,log10(abs(FT_fft_CI).^2));
xlabel('Range')
% ylim([0 1.8])
ylabel('Normalized Power (u.W) ')
title('Coherent Integration of every Chirp')
Index_rangemax = find(abs(FT_fft_CI)==max(abs(FT_fft_CI)));
text(axis_range(Index_rangemax),abs(FT_fft_CI(Index_rangemax)).^2, [' Range measured: ' num2str(axis_range(Index_rangemax)) ,' (',num2str(R),')'] );

%% Phase Extract of signal before shift
I = real(FT_fft(:,Index_rangemax));
Q = imag(FT_fft(:,Index_rangemax));
Phi = atan(I./Q);
Phi = unwrap_def(Phi,pi/2);
Trajectory = Phi*lambda/(4*pi);
Trajectory = detrend(Trajectory,0);

% figure;
% plot(Phi)
figure;
plot(t_slow ,Trajectory);
xlabel('Time (s)')
ylabel('Amplitude (m)')
title('Trajectory of Mition Before Shift')
% ylim([-1.5*Am 1.5*Am])

Amp_measure = max(Trajectory) - min(Trajectory);

%% Frequency Shift 
ShiftStepRes = 4;
ShiftV = -ShiftStepRes/2:ShiftStepRes/2;
f_shiftstep =1/ShiftStepRes * delta_f;
% f_shiftstep = 200;

Peaks = zeros(1,ShiftStepRes+1);
C = zeros(1,ShiftStepRes+1);
H = zeros(1,ShiftStepRes+1);



for i = -ShiftStepRes/2:ShiftStepRes/2
    X_fshift = fft(ifft(FT_fft_CI).*exp(1i*2*pi*f_shiftstep*t_fast*i));
    X_Ashift = abs(X_fshift);

    figure;
    plot(axis_range,X_Ashift,'LineWidth',1.5);
    xlabel('Fre') 
    TitleName = ['Coherent Integration of every Chirp After Shift ',num2str(i/ShiftStepRes*delta_range)];
    title(TitleName,'FontSize',12)
    ylim([0 0.035])

    Index_max = find(abs(X_fshift)==max(abs(X_fshift)));
    text(axis_range(Index_max), abs(X_fshift(Index_max)), ['(',num2str(axis_range(Index_max)), ',',num2str(abs(X_fshift(Index_max))) ,')'],'FontSize',12);
    
    X_Afactor = X_Ashift/max(X_Ashift)*10;
    C(i+ShiftStepRes/2+1) = var(X_Afactor)/mean(X_Afactor);
    H(i+ShiftStepRes/2+1) = log(sum(X_Afactor)) - 1/sum(X_Afactor) * sum( X_Afactor.*log(X_Afactor));
    Peaks(i+ShiftStepRes/2+1) = max(X_Ashift);

end

Index_Cmax = find(C == max(C));
figure;
plot(C)
title("Quality factor C")
text(Index_Cmax,C(Index_Cmax),['Max Index is',num2str(Index_Cmax)])

Index_Hmin = find(H == min(H));
figure;
plot(H)
title("Quality factor H")
text(Index_Hmin,H(Index_Hmin),[' Min Index is ',num2str(Index_Hmin)])

Index_Peaks = find(Peaks == max(Peaks));
figure;
plot(Peaks,'LineWidth',1.5)
title("Peak Search",'FontSize',12)

text(Index_Peaks,Peaks(Index_Peaks),[' Max Index is ',num2str(Index_Peaks)],'FontSize',12)
%% Add Phase Ramp With Signal

ChirpArray_Shift = ChirpArray.*exp(1i*2*pi*f_shiftstep*t_fast*ShiftV(Index_Peaks));
Range_fft =fft(ChirpArray_Shift')';
figure;
mesh(axis_range,t_slow,abs(Range_fft))
axis tight;
% view(0,90);
colorbar;
xlabel('Range')
ylabel('Slow Time Index')
title('Range-FFT')

%%  CPI
FT_fft_CI = sum((Range_fft))/Ns_chirp/Nc;

figure;
plot(axis_range,log10(abs(FT_fft_CI).^2));
xlabel('Range')
% ylim([0 1.8])
ylabel('Normalized Power (u.W) ')
title('Coherent Integration of every Chirp')
Index_rangemax = find(abs(FT_fft_CI)==max(abs(FT_fft_CI)));
text(axis_range(Index_rangemax),abs(FT_fft_CI(Index_rangemax)).^2, [' Range measured: ' num2str(axis_range(Index_rangemax)) ,' (',num2str(R),')'] );

%% Phase extraction of signal after shift
I = real(Range_fft(:,Index_rangemax));
Q = imag(Range_fft(:,Index_rangemax));
Phi = atan(I./Q);
Phi = unwrap_def(Phi,pi/2);
Trajectory = Phi*lambda/(4*pi);
Trajectory = detrend(Trajectory,0);

% figure;
% plot(Phi)
figure;
plot(t_slow ,Trajectory);
xlabel('Time (s)')
ylabel('Amplitude (m)')
title('Trajectory of Mition After Shift')
% ylim([-1.5*Am 1.5*Am])


%% constellation chart
figure;
plot(I,Q);
axis equal



