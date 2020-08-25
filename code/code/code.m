clear
clc
%% iMPULSE

sr = 44100;

durata_l = 1;
durata_r = 1;
lunghezza = 1 * 44100;

T_sign_l = linspace(0,durata_l, lunghezza);
Y_l = zeros(1, lunghezza);
Y_l(lunghezza/2) = 1;
% Y_l(lunghezza/2) = 1;

T_sign_r = linspace(0, durata_r, lunghezza);
Y_r = zeros(1, lunghezza);
% Y_r(lunghezza/2) = 1;
Y_r(lunghezza/2 + 882) = 1;

%% CHORD TEST

sr = 44100;

nfft = 4096;
fin_f = hamming(nfft);
stepSize = nfft/4;
overlap_f = nfft - stepSize;

durata_l = 3;
durata_r = 3;
lunghezza = 3 * sr;

f0 = 261.6;
f2 = 329.6;
f3 = 392;

f0r = 260.6;
f2r = 330.4;
f3r = 393;

Y_l = zeros(1, lunghezza);
Y_r = zeros(1, lunghezza);
T_sign_l = linspace(0,durata_l, lunghezza);
T_sign_r = linspace(0, durata_r, lunghezza);

for i = 2048 : lunghezza
    if i < sr
        Y_l(i) = 0.5*sin(2*pi*f0*(i-1)/sr) + 0.3*sin(2*pi*2*f0*(i-1)/sr) + 0.15*sin(2*pi*3*f0*(i-1)/sr) + 0.08*sin(2*pi*4*f0*(i-1)/sr);
    end
    if i < sr*2 && i >= sr
        Y_l(i) = 0.7*sin(2*pi*f2*(i-1)/sr) + 0.3*sin(2*pi*2*f2*(i-1)/sr) + 0.15*sin(2*pi*3*f2*(i-1)/sr) + 0.08*sin(2*pi*4*f2*(i-1)/sr);
    end
    
    if i < sr*3 && i >= sr*2
        Y_l(i) = 0.8*sin(2*pi*f3*(i-1)/sr) + 0.3*sin(2*pi*2*f3*(i-1)/sr) + 0.15*sin(2*pi*3*f3*(i-1)/sr) + 0.08*sin(2*pi*4*f3*(i-1)/sr);
    end
end
for i = 2048+1024 : lunghezza
    if i < sr + 1024
        Y_r(i) = 0.5*sin(2*pi*f0r*(i-1)/sr) + 0.3*sin(2*pi*2*f0r*(i-1)/sr) + 0.15*sin(2*pi*3*f0r*(i-1)/sr) + 0.08*sin(2*pi*4*f0r*(i-1)/sr);
    end
    if i < sr*2 + 1024 && i >= sr + 1024
        Y_r(i) = 0.7*sin(2*pi*f2r*(i-1)/sr) + 0.3*sin(2*pi*2*f2r*(i-1)/sr) + 0.15*sin(2*pi*3*f2r*(i-1)/sr) + 0.08*sin(2*pi*4*f2r*(i-1)/sr);
    end
    
    if i < sr*3 && i >= sr*2 + 1024
        Y_r(i) = 0.8*sin(2*pi*f3r*(i-1)/sr) + 0.3*sin(2*pi*2*f3r*(i-1)/sr) + 0.15*sin(2*pi*3*f3r*(i-1)/sr) + 0.08*sin(2*pi*4*f3r*(i-1)/sr);
    end
end
    
%% SIN TEST
clc

% SIN TEST
sr = 44100;
durata_l = 2;
lung = sr*2;
% f = 500*sr/8192;
f = 880;
f1 = 880.5;
T_sign_l = linspace(0,durata_l,lung); 

for i = 1 : lung
   Y_l(i) = 0.5*cos(2*pi*f*(i - 1)/sr); 
   Y_r(i) = 0.5*cos(2*pi*f1*(i - 1)/sr);
end
Export(:,1) = Y_l;
Export(:,2) = Y_r;
audiowrite('Output.wav', Export, sr);

%% CHIRP TEST
clc

% CHIRP TEST
sr = 44100;
durata_l = 3;
lung = durata_l * sr;
T_sign_l = linspace(0,durata_l,lung);
nfft = 2048;
f0 = 200;
f1 = 500;
Y_l = chirp(T_sign_l, f0, durata_l -1 , f1);

audiowrite('Output.wav', Y_l, sr);

%% IMPORT
clc
% 
[Y_L, ~] = audioread('ch.wav');             %  Carico primo file
[Y_R, sr] = audioread('ch.wav');            %  Carico secondo file

%   CONTROLLO MONO
if size(Y_L,2) > 1 
    Y_l = Y_L(:,1);
else
    Y_l = Y_L;
end
if size(Y_R,2) > 1
    Y_r = Y_R(:,2);
else
    Y_r = Y_R;
end
% 
Y_l = Y_l(2048:end-11025);
Y_r = Y_r(2048:end-11025);
% 
% start = 0.80;
% interval = 1.1;
% ending = floor(start*sr) + floor(interval*sr);
% step = floor(interval*sr);
% Y_l = Y_l(floor(start*sr):ending);
% Y_r = Y_r(floor(start*sr):ending);

%   LUNGHEZZA SEGNALE LEFT
lunghezza_l = size(Y_l,1);                          %   Lunghezza in campioni (samples)
durata_l = lunghezza_l/sr;                          %   Durata in secondi (s)

T_sign_l = linspace(0,durata_l,lunghezza_l);        %   Creo vettore temporale

%   LUNGHEZZA SEGNALE RIGHT
lunghezza_r = size(Y_r,1);                          %   Lunghezza in campioni (samples)
durata_r = lunghezza_r/sr;                          %   Durata in secondi (s)

T_sign_r = linspace(0,durata_r,lunghezza_r);        %   Vettore temporale

%   FINESTRATURE DOMINIO DEL TEMPO
finestra_t = 4096;                                   %   Time Domain Window
overlap_t = finestra_t - (finestra_t/4);                 %   Overlap

%   FINESTRATURE DOMINIO FREQUENZA

nfft = 1024;                                %   Freq Domain Window == Finestra analisi
stepSize = nfft/4;
overlap_f =  nfft - stepSize;                 %   Overlap
fin_f = hamming(nfft);                        %   Window functio

%%   ANALISI NEL DOMINIO DEL TEMPO ( RMS - PICCO - CREST FACTOR)
clc

%   DIVIDO IL SEGNALE IN FRAME
Y_left = buffer(Y_l,finestra_t,overlap_t);
Y_right = buffer(Y_r,finestra_t,overlap_t);

%   RMS
rms_left = mag2db(rms(Y_left));
rms_right = mag2db(rms(Y_right));
rms_maxLeft = max(rms_left);
rms_maxRight = max(rms_right);
rms_minLeft = min(rms_left);
rms_minRight = min(rms_right);
rms_meanLeft = mean(rms_left);
rms_meanRight = mean(rms_right);
rms_stdLeft = std(rms_left, 0);
rms_stdRight = std(rms_right, 0);

rms_diff = abs(rms_left - rms_right);
% % PICCO
picco_left = mag2db(picco(Y_left));
picco_right = mag2db(picco(Y_right));
picco_maxLeft = max(picco_left);
picco_maxRight = max(picco_right);
picco_minLeft = min(picco_left);
picco_minRight = min(picco_right);
picco_meanLeft = mean(picco_left);
picco_meanRight = mean(picco_right);
picco_stdLeft = std(picco_left, 0);
picco_stdRight = std(picco_right, 0);

picco_diff = abs(picco_left - picco_right);
% % CREST FACTOR
crestF_left = (crestFact(picco_left, rms_left));
crestF_right = (crestFact(picco_right, rms_right));
crestF_maxLeft = max(crestF_left);
crestF_maxRight = max(crestF_right);
crestF_minLeft = min(crestF_left);
crestF_minRight = min(crestF_right);
crestF_meanLeft = mean(crestF_left);
crestF_meanRight = mean(crestF_right);
crestF_stdLeft = std(crestF_left, 0);
crestF_stdRight = std(crestF_right, 0);

crest_diff = abs(crestF_left - crestF_right);
 
T_din_l = linspace(0,durata_l,size(Y_left,2));          %   Vettore temporale
T_din_r = linspace(0,durata_r,size(Y_right,2));


%% 

subplot(2,1,1);
plot(T_din_l, (rms_diff), 'linewidth', 1.5);
title('RMS difference');
xlabel('Tempo (s)');
ylabel('Amp (dB)');
subplot(2,1,2);
plot(T_din_l, (picco_diff), 'linewidth', 1.5);
title('Picco difference');
xlabel('Tempo (s)');
ylabel('Amp (dB)');
axis([0 1.8 5.0005 5.001]);
% subplot(3,1,3);
% plot(T_din_l, (crest_diff), 'linewidth', 1);
% title('Crest Factor difference');
% xlabel('Tempo (s)');
% ylabel('Amp (dB)');

diff_mean_picco = mean(picco_diff(1:end));
diff_std_picco = std(picco_diff(1:end));

diff_mean_rms = mean(rms_diff(1:end));
diff_std_rms = std(rms_diff(1:end));

diff_mean_crest = mean(crest_diff(1:end));
diff_std_crest = std(crest_diff(1:end));

%% MS CODING
clc

[Mid, Side] = ms(Y_l, Y_r);

%% TEST MS
% clc
% 
% figure;
% subplot(1,2,1)
%     plot( Y_l(44100:88200), Y_r(44100:88200), '.', 'markersize', 1 )
%     line([0 0],[-1 1],'linestyle',':','color','k');
%     line([-1 1],[0 0],'linestyle',':','color','k');
%     grid on; 
%     title('Left vs Right');
%     xlabel('Dati canale L');
%     ylabel('Dati canale R');
%     axis([-1, 1, -1, 1]);
% subplot(1,2,2)
%     plot( Side(44100:88200), Mid(44100:88200), '.', 'markersize',1);
%     line([-1 1],[-1 1],'linestyle',':','color','k');
%     line([-1 1],[1 -1],'linestyle',':','color','k');
%     grid on; 
%     title('Mid vs Side');
%     xlabel('Dati canale S');
%     ylabel('Dati canale M');
%     axis([-1, 1, -1, 1]);

%% ANALISI NEL DOMINIO DELLA FREQUENZA : STFT 
clc

% SEGNALE LEFT
[STFT_l, F_l, T_sp_l] = spectrogram(Y_l, fin_f, overlap_f, nfft, sr);
STFTmag_l = abs(STFT_l)./(nfft).*2;
STFTfase_l = angle(STFT_l);

%SEGNALE RIGHT
[STFT_r, F_r, T_sp_r] = spectrogram(Y_r, fin_f, overlap_f, nfft, sr);
STFTmag_r = abs(STFT_r)./(nfft).*2;
STFTfase_r = angle(STFT_r);

% [Y_left2, time] = istft(STFT_l, nfft, stepSize, nfft, sr);

% %% LAG TIME
% clc
% % 
% % Y_left2 = buffer(Y_l, nfft, overlap_f);
% % Y_left2 = buffer(Y_l, nfft, overlap_f);
% 
% 
% nfftcorr = 4096;
% fin_fcorr = hamming(nfftcorr);
% stepSizecorr = nfft/4;
% overlap_fcorr = nfftcorr - stepSize;
% 
% % SEGNALE LEFT
% [STFT_lcorr, F_lcorr, T_sp_lcorr] = spectrogram(Y_l, fin_fcorr, overlap_fcorr, nfftcorr, sr);
% STFTmag_lcorr = abs(STFT_lcorr)./(nfftcorr);
% % STFTfase_l = angle(STFT_l);
% 
% %SEGNALE RIGHT
% [STFT_rcorr, F_rcorr, T_sp_rrcorr] = spectrogram(Y_r, fin_fcorr, overlap_fcorr, nfftcorr, sr);
% STFTmag_rcorr = abs(STFT_rcorr)./(nfftcorr);
% % STFTfase_r = angle(STFT_r);
% 
% 
% for c = 1 : length(T_sp_lcorr)
%     [Y_left2, ~] = istft(STFT_lcorr(:,c),nfftcorr, stepSizecorr, nfftcorr, sr);
%     [Y_right2, ~] = istft(STFT_rcorr(:,c),nfftcorr, stepSizecorr, nfftcorr, sr);
% %     Y_left2 = Y_left2./fin_fcorr;
% %     Y_right2 = Y_right2./fin_fcorr;
% %     Y_left2 = rms(Y_left2);
% %     Y_right2 = rms(Y_right2);
%     [corr, lag] = xcorr(Y_left2, Y_right2);
%     [~, I] = max(abs(corr));
%     lagD(c) = lag(I)/sr;
% end


%% Frequency Estimation
clc

FE_l = diffPhase(STFTfase_l, nfft, nfft/stepSize, sr);
FE_r = diffPhase(STFTfase_r, nfft, nfft/stepSize, sr);

%% PICCHI
clc

thresh = -42.7 + 8;

Picchi_l = Peaks(STFTmag_l, thresh, 1);
Picchi_r = Peaks(STFTmag_r, thresh, 1);

%% TEST PICCHI
clc

figure;
semilogx(F_l,mag2db(STFTmag_r(:,2)));
hold on;
% semilogx(F_l,massi(:,200));
hold on;
semilogx(F_l,Picchi_r(:,2),'x');

%% MASSIMI
clc

%   PICCHI MASSIMI
[picchiMassimi_l, valoreMassimi_l] = estraiMax(Picchi_l, STFTmag_l, 1);
[picchiMassimi_r, valoreMassimi_r] = estraiMax(Picchi_r, STFTmag_r, 1);

%% TEST MASSIMI
clc

figure;
semilogx(F_l,mag2db(STFTmag_r(:,10)), 'linewidth', 1);
hold on;
% semilogx(F_l,Mob_l(:,200));
hold on;
semilogx(F_l,valoreMassimi_r(:,10),'*k');
hold off;
xlabel('Frequenza (Hz)');
ylabel('Ampiezza(dB)');
legend('FFT','Picco');
grid on

%% DIFFERENZA DI FASE E UNWRAPPING
clc

%   CALCOLO DIFERENZA DI FASE E PHASE UNWRAPPING L
wrap_l = zeros(size(STFTfase_l,1),size(STFTfase_l,2));
for c = 1: size(STFTfase_l,1)
    wrap_l(c,:) = unwrap(STFTfase_l(c,:));
end

%   CALCOLO DIFERENZA DI FASE E PHASE UNWRAPPING R
wrap_r = zeros(size(STFTfase_r,1),size(STFTfase_r,2));
for c = 1: size(STFTfase_r,1)
    wrap_r(c,:) = unwrap(STFTfase_r(c,:));
end

diff_l = diff(wrap_l,1,2);
diff_r = diff(wrap_r,1,2);

diffW_l = unwrap(wrap_l,[],1);
diffW_r = unwrap(wrap_r,[],1);

%% INTERPOLAZIONE PARABOLICA
clc

%   INTERPOLAZIONE PARABOLICA
[locPicchi_l, valPicchi_l] = interParab(picchiMassimi_l, STFTmag_l);
valPicchi_l = mag2db(valPicchi_l);
[locPicchi_r, valPicchi_r] = interParab(picchiMassimi_r, STFTmag_r);
valPicchi_r = mag2db(valPicchi_r);

%% F0 Detection
clc

% Cepstrum Based - DETECTED

% SEGNALE LEFT
CpCoeff_l = Cepstrum(STFT_l);
% CpCoeff_l = CpCoeff_l(1:end,:);
F0_l = PitchCepstrum(CpCoeff_l, sr);
F0_l = F0_l./2;
F0_l = medfilt1(F0_l, 20);
%SEGNALE RIGHT
CpCoeff_r = Cepstrum(STFT_r);
CpCoeff_r = CpCoeff_r(2:end,:);
F0_r = PitchCepstrum(CpCoeff_r, sr);
F0_r = F0_r./2;
F0_r = medfilt1(F0_r, 20);

figure('Name','F0 Estimation with Cepstrum','NumberTitle','off');
% imagesc(T_sp_l, F_l(1:nfft/8),mag2db(STFTmag_r(1:nfft/8,:)));
% c = colorbar; c.Label.String = 'Magnitude dB';
% axis xy;
% hold on
% semilogy(T_sp_l,F0_r);
% hold on;
semilogy(T_sp_l,F0_l,'linewidth', 2);
hold on;
semilogy(T_sp_r,F0_r,'linewidth', 2);
grid on
legend Left Right
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
% print('-dpdf', '01 - f0_estimate');

%% FUNDAMENTAL DIFFERENCE

F0diff = (F0_l ./ F0_r);
F0diff = abs(log(F0diff)./log(2.^(1/1200)));

F0_lmedian = median(F0_l);
F0_rmedian = median(F0_r);

F0_lrmedian = median(F0diff);

% 
figure;
plot(T_sp_l, F0diff, 'linewidth', 1);
xlabel('Tempo (s)');
ylabel('Cents');

% 
F0mean = mean(F0diff);
F0std = std(F0diff, 0);
%% SINE TRACKING PHASE DIFFERENCE
clc

%SEGNALE LEFT
picchiMassimi_l = NaN2zero(picchiMassimi_l);
[Out_l] = sinTracking(sr, nfft, picchiMassimi_l, FE_l, 15);
% [Out_l] = deleteEl(Out_l, F0_l);
[Out_l] = cleanTrack(Out_l, 10);
Out_l = zero2NaN(Out_l);

%SEGNALE RIGHT

picchiMassimi_r = NaN2zero(picchiMassimi_r);
[Out_r] = sinTracking(sr, nfft, picchiMassimi_r, FE_r, 15);
% [Out_r] = deleteEl(Out_r, F0_r);
[Out_r] = cleanTrack(Out_r, 10);
Out_r = zero2NaN(Out_r);

%% SIN TRAKING PARABOLIC INTERPOLATION
clc

% SEGNALE LEFT
locPicchi_l = NaN2zero(locPicchi_l);
[outParab_l] = sinTrackingParab(sr, nfft, locPicchi_l, 30);
% [outParab_l] = deleteEl(outParab_l, F0_l);
[outParab_l] = cleanTrack(outParab_l, 5);
outParab_l = zero2NaN(outParab_l);

%SEGNALE RIGHT
locPicchi_r = NaN2zero(locPicchi_r);
[outParab_r] = sinTrackingParab(sr, nfft, locPicchi_r, 30);
% [outParab_r] = deleteEl(outParab_r, F0_r);
[outParab_r] = cleanTrack(outParab_r, 5);
outParab_r = zero2NaN(outParab_r);


%% TEST SINE TRAKING
clc
%SEGNALE LEFT

figure('Name','Tracking Fase','NumberTitle','off');
imagesc(T_sp_l, F_l(1:floor(nfft/15)),mag2db(STFTmag_l(1:floor(nfft/15),:)));
% c = colorbar; c.Label.String = 'Magnitude dB';
axis xy;
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
hold on;
semilogy(T_sp_l,Out_r, 'black','LineWidth',1);
hold off;
xlabel('Tempo (s)');
ylabel('Frequenza(Hz)');
grid on
% 
% figure('Name','Traking Parabola','NumberTitle','off');
% imagesc(T_sp_l, F_l(10:floor(nfft/40)),mag2db(STFTmag_l(10:floor(nfft/40),:)));
% c = colorbar; c.Label.String = 'Magnitude dB';
% % colormap('jet');
% axis xy;
% hold on;
% plot(T_sp_l,outParab_l, 'black','LineWidth',1);
% hold off
% xlabel('Tempo (s)');
% ylabel('Frequenza(Hz)');

%SEGNALE RIGHT

% figure('Name','Traking Fase R','NumberTitle','off');
% imagesc(T_sp_r,  F_r(1:nfft/8),mag2db(STFTmag_r(1:nfft/8,:)));
% c = colorbar; c.Label.String = 'Magnitude dB';
% axis xy;
% hold on;
% plot(T_sp_r,Out_r, 'black','LineWidth',1);
% hold off;
% 
% figure('Name','Traking Parabola','NumberTitle','off');
% imagesc(T_sp_r, F_r(1:nfft/8,:),mag2db(STFTmag_r(1:nfft/8,:)));
% c = colorbar; c.Label.String = 'Magnitude dB';
% axis xy;
% hold on;
% plot(T_sp_r,outParab_r, 'black','LineWidth',1);
% hold off

%% OTHER TEST

figure('Name','harmonic tracking');
semilogy(T_sp_l,Out_l, 'or');
hold on;
semilogy(T_sp_r,Out_r, '-k', 'linewidth', 2);
hold off;
grid on
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');

%%

% minimo = +inf;
% cent = log(2.^(1/1200));
% threshcent = 400;
% 
% Sinistro = zeros(size(Out_l,1), length(T_sp_l));
% Destro = zeros(size(Out_r,1), length(T_sp_r));
% match = NaN(1, length(T_sp_l));
% 
% b = 1;
% 
% for i =1 : T_sp_l
%    uno = 1;
%    due = 1;
%    for j =1 : size(Out_l,1)
%        loc = 0;
%        for k = 1 : size(Out_r,1)
%        match(k) = abs(log(Out_l(j,i) ./ Out_r(k,i))./cent);
%             if match(k) < minimo &&  match(k) < threshcent
%                 minimo = match(k);
%                 loc = k;
%             else
%                 continue;
%             end
%        end
%         if loc == 0
%             continue;
%         end
%     Sinistro(uno,i) = Out_l(j,i);
%     Destro(due,i) = Out_r(loc,i);
%     Out_r(loc,i) = inf;
%     uno = uno + 1;
%     due = due + 1;
%     minimo = +inf;
%    end
% end

%% HARMONIC DIFFERENCE
clc
nharm = 1;
cent = log(2.^(1/1200));

% Phase Difference Method
for harm = 1 : nharm
    PHdiffHarm(harm,:) = abs(log(Out_l(harm,:) ./ Out_r(harm,:)) ./ cent);
    frequenza(harm,:) = abs(Out_l(harm,:) - Out_r(harm,:));
end

PHminDiffHarm = min(PHdiffHarm, [], 2);
PHmaxDiffHarm = max(PHdiffHarm, [], 2);
PHdiffHarm = NaN2zero(PHdiffHarm);
PHmedianDiffHarm = mean(PHdiffHarm, 2);
PHstdDiffHarm = std(PHdiffHarm, 0, 2);
PHdiffHarm = zero2NaN(PHdiffHarm);

ffre = mean(frequenza(2:end));

PHmm = mean(PHmedianDiffHarm);
PHstd = std(PHmedianDiffHarm);

%%
subplot(2,1,1);
plot(T_sp_l, PHdiffHarm, 'linewidth', 2);
xlabel('Tempo (s)');
ylabel('Cents');
title('Cent difference');
subplot(2,1,2);
plot(T_sp_l, frequenza, 'linewidth', 2);
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
title('Frequency difference');

%% PERCUSSIVE FEATURE DETECION
clc

nfftPFD = 512;
fin_fPFD = hann(nfftPFD);
stepSizePFD = nfftPFD/8;
overlap_fPFD = nfftPFD - stepSizePFD;

% SEGNALE LEFT
[STFT_lPFD, F_lPFD, T_sp_lPFD] = spectrogram(Y_l, fin_fPFD, overlap_fPFD, nfftPFD, sr);
STFTmag_lPFD = abs(STFT_lPFD)./(nfftPFD);

%SEGNALE RIGHT
[STFT_rPFD, F_rPFD, T_sp_rPFD] = spectrogram(Y_r, fin_fPFD, overlap_fPFD, nfftPFD, sr);
STFTmag_rPFD = abs(STFT_rPFD)./(nfftPFD);
% 
% [YPFDL, ~] = istft(STFT_lPFD(:,2),nfftPFD, stepSizePFD, nfftPFD, sr);
% [YPFDR, ~] = istft(STFT_rPFD(:,2),nfftPFD, stepSizePFD, nfftPFD, sr);

onset_l = PFD(STFTmag_lPFD, 0.2);
onset_r = PFD(STFTmag_rPFD, 0.2);

% plot(T_sp_lPFD(1:end-1),onset_r)

onPeak_l = Peaks(onset_l, 11.5, 0);
onPeak_r = Peaks(onset_r, 11.5, 0);

[onPos_l, onMax_l]  = estraiMax(onPeak_l, onset_l, 0);
[onPos_r, onMax_r]  = estraiMax(onPeak_r, onset_r, 0);

onTime_l = linspace(0,durata_l,length(onMax_l));
onTime_r = linspace(0,durata_r,length(onMax_r));

% LEFT 1° riga, RIGHT 2° riga

u = 1;
i = 1;

for c = 1 : length(onMax_l)
    if isnan(onMax_l(c)) && isnan(onMax_r(c))
        continue;
    end
    if ~isnan(onMax_l(c))
        newL(1,i) = max(1, c - 1);
        i = i + 1;
    end
    if ~isnan(onMax_r(c))
        newL(2,u) = max(1, c - 1);
        u = u + 1;
    end
end

figure('Name', 'Plot Correlaztion');
plot(onTime_l, onPos_l + 1, 'o');
hold on
plot(onTime_l, onPos_r + 1, 'xk');
hold off
legend Left Right
grid on
xlabel('Tempo (s)');
print('-dpdf', '02 - onset.pdf');

figure;
plot(onTime_l,onset_r,'linewidth', 1)
hold on;
plot(onTime_l,onset_l, 'linewidth', 1)
legend Destro Sinistro
title('Percussive Feature Detection');
xlabel('Tempo (s)');
ylabel('bin');

% figure('Name', 'Plot Onset');
% hold on
% plot(onTime_l, onMax_l, 'o');
% hold on 
% plot(onTime_r, onMax_r, 'x');
% hold on 
% plot(onTime_l, onPeak_l);
% hold on;
% plot(onTime_r, onPeak_r);
% hold off
% legend Left Right
% grid on

%%
figure;
subplot(2,1,1);
plot(T_sign_l, Y_l,'black','linewidth', 1);
title('Waveform');
xlabel('Tempo (s)');
ylabel('Ampiezza');
grid on
subplot(2,1,2);
plot(onTime_l,onset_l,'black','linewidth', 1);
title('Percussive Feature Detection');
xlabel('Tempo (s)');
ylabel('bin');
grid on;

%% LAG TIME ONSET 

nframe = 1;
newL = newL(:, 1: nframe);

for c = 1 : size(newL,2)
%     L = Y_l(newL(1, c)*stepSizePFD : min(newL(1, c)*stepSizePFD + nfftPFD));
%     R = Y_r(newL(2, c)*stepSizePFD : min(newL(2, c)*stepSizePFD + nfftPFD));
%     [crossLag, lags] = xcorr(L, R);
%     [~, I] =Per max(abs(crossLag));
%     lagD(c) = lags(I)/sr;
%     lagD(c) = lagD(c) - (newL(1,c) - newL(2, c)) * nfftPFD/sr;
    lagD(c) = T_sp_lPFD(newL(1, c)) - T_sp_rPFD(newL(2, c));
    lagDabs(c) = abs(T_sp_lPFD(newL(1, c)) - T_sp_rPFD(newL(2, c)));
end

%%

lagmean = mean(lagDabs);
lagstd = std(lagDabs);


%% BEST LPC ORDER

yl = zeros(nfft, length(T_sp_l));
yr = zeros(nfft, length(T_sp_l));

for c = 1 : length(T_sp_l)
    [yl(:,c), ~] = istft(STFT_l(:,c), nfft, stepSize, nfft, sr);  
    [yr(:,c), ~] = istft(STFT_r(:,c), nfft, stepSize, nfft, sr); 
end

ncoef = 1 + round(sr / 1000);
Order_l = Ordine(yl, ncoef);
Order_r = Ordine(yr, ncoef);

Or = floor(0.5*(Order_l + Order_r)/2);

if Or >= 23
    Or = 23;
end

%% LPC Envelope
clc
pthres = 0.005;
LpC_l = zeros(Or + 1, length(T_sp_l));
LpC_r = zeros(Or + 1, length(T_sp_l));
Freq_l = NaN(floor(Or), length(T_sp_l));
Freq_r = NaN(floor(Or), length(T_sp_l));
 
for c = 1 : length(T_sp_l)  
%     [a_l(:,c), residual_l, gain_l] = LinearPredCod(yl(:,c), Order_l);
    [LpC_l(:,c), e_l] = aryule(yl(:,c), Or);
    r_l = roots(LpC_l(:,c));
    ramp_l = real(r_l);
    r_l(imag(r_l) < pthres) = nan;
    r_l = angle(r_l);
    r_l = sort(r_l.*(sr/(2*pi)));
    Freq_l(:,c) = r_l;
%     [a_r(:,c), residual_r, gain_r] = LinearPredCod(yr(:,c), Order_r);
    [LpC_r(:,c), e_r] = aryule(yr(:,c), Or);
    r_r = roots(LpC_r(:,c));
    ramp_r = real(r_r);
    r_r(imag(r_r) < pthres) = nan;
    r_r = angle(r_r);
    r_r = sort(r_r.*(sr/(2*pi)));
    Freq_r(:,c) = r_r;  
    [H_l(:,c), P_l] = freqz(sqrt(e_l), LpC_l(:,c), length(F_l), sr);
    [H_r(:,c), P_r] = freqz(sqrt(e_r), LpC_r(:,c), length(F_r), sr);
    H_l(:,c) = abs(H_l(:,c));
    H_r(:,c) = abs(H_r(:,c));
    LSD(:,c) = (sqrt(mean((mag2db(H_l(:,c)) - mag2db(H_r(:,c))).^2)));
end

LSDmean = mean(LSD);

%% FIRST N POLES MATCHES

minimo = +inf;
numpoles = 6;
cent = log(2.^(1/1200));
threshcent = 800;

polesLeft = zeros(numpoles, length(T_sp_l));
polesRight = zeros(numpoles, length(T_sp_r));
mismatch = NaN(1, floor(Or/2));

b = 1;

for i = 1 : length(T_sp_l)
        uno = 1;
        due = 1;
    for j = 1 : floor(Or/2)
        loc = 0;
        for k = 1 : floor(Or/2)
            mismatch(k) = abs(log(Freq_l(j,i) ./ Freq_r(k,i))./cent);
            if mismatch(k) < minimo &&  mismatch(k) < threshcent
                minimo = mismatch(k);
                loc = k;
            else
                continue;
            end
        end
        if loc == 0
            continue;
        end
        polesLeft(uno,i) = Freq_l(j,i);
        polesRight(due,i) = Freq_r(loc,i);
        Freq_r(loc,i) = inf;
        uno = uno + 1;
        due = due + 1;
        minimo = +inf;
    end
end

%% MEDIA E STD DELLE FORMANTI

% Freq_l = NaN2zero(Freq_l);

media_f = mean(polesLeft(1:floor(size(polesLeft, 1)/2),:), 2);
devstd = std(Freq_l(1:floor(size(polesRight, 1)/2),:), 0, 2);

media_l = mean(H_l, 2);
devstd_l = std(H_l, 0, 2);
media_r = mean(H_r, 2);
devstd_r = std(H_r, 0, 2);

mediapole_l = mean(polesLeft(1:numpoles,:), 2);
devstdpole_l = std(polesLeft(1:numpoles,:), 0, 2);
mediapole_r = mean(polesRight(1:numpoles,:), 2);
devstdpole_r = std(polesRight(1:numpoles,:), 0, 2);

figure('Name', 'Media Formanti');
semilogx(P_l, media_l - 6000);
hold on
semilogx(P_r, media_r - 6000);
hold on
semilogx(mediapole_l, ones(length(mediapole_l),1), 'or');
hold on
semilogx(mediapole_r, ones(length(mediapole_r),1), '*b');
figure('Name', 'Deviazione Standard  Formanti');
semilogx(P_l, devstd_l);
hold on
semilogx(P_r, devstd_r);

%mediapole
%%

figure('Name','Tracking Fase','NumberTitle','off');
imagesc(T_sp_l, F_l(1:floor(nfft/2)),mag2db(STFTmag_l(1:floor(nfft/2),:)));
colormap('jet');
axis xy;
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
hold on;
semilogx(T_sp_l, polesLeft(1:numpoles,:), '*k');
%% TEST FORMANT LPC BASED
clc
figure;
semilogx(P_l, mag2db(abs(H_l(:,2))),'r','Linewidth', 1.5);
hold on;
semilogx(P_r, mag2db(abs(H_r(:,2))),'k','Linewidth', 1.5);
hold on;
semilogx(F_l, mag2db(STFTmag_l(:,2)),'r');
hold on;
semilogx(F_r, mag2db(STFTmag_r(:,2)),'k');
hold on
semilogx(polesLeft(1:numpoles,2), ones(numpoles,1), 'or');
hold on
semilogx(polesRight(1:numpoles,2), ones(numpoles,1), '*k');
hold off;
legend ('Lpc Left','Lpc Right', 'FFT Left','FFT Right', 'Formanti Left','Formanti Right');
grid on
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
% print('-dpdf', 'formant.pdf');

%%

cent = log(2.^(1/1200));

for formant = 1 : numpoles
    diffPole(formant,:) = abs(log(polesLeft(formant,:) ./ polesRight(formant,:)) ./ cent);
end

diffmeanPole = mean(diffPole,2);
diffstdPole = std(diffPole,0,2);

ffmm = mean(diffmeanPole);
ffstd = std(diffmeanPole);
%% TEST INTERPOLAZIONE PARABOLICA
clc

figure;
semilogx(F_l,mag2db(STFTmag_l(:,5)),'linewidth',1);
hold on;
semilogx(((sr/2)*(locPicchi_l(:,5) - 1)/length(F_l)),valPicchi_l(:,5),'xk');
hold off; 
xlabel('Frequenza(Hz)');
ylabel('Ampiezza(dB)');
legend('FFT', 'picco');
grid on

%% TEST Differenza fase
clc

figure;
semilogx(F_l,mag2db(STFTmag_l(:,5)),'linewidth',1);
hold on;
semilogx(Out_l(:,5),zeros(1,length(Out_l)) - 19,'xk');
hold off; 
xlabel('Frequenza(Hz)');
ylabel('Ampiezza(dB)');
legend('FFT', 'picco');
grid on

%%

plot(T_sp_l,Out_l,'LineWidth',2);
xlabel('Tempo (s)');
ylabel('Frequenza (Hz)');
hold on
plot(T_sp_l,outParab_l,'LineWidth',2);
grid on
legend ('Differenza di fase','Interpolazione parabolica');

%% PLOT SPETTROGRAMMI L & R

%Forma d'onda e Spettrogramma L
figure('Name','Forma d^onda e Spettrogramma','NumberTitle','off');
subplot(2,1,1);
plot(T_sign_l, Y_l);
xlabel('Tempo (s)');
ylabel('Ampiezza');
grid on
subplot(2,1,2);
imagesc(T_sp_l,F_l(1:floor(end/5)),mag2db(STFTmag_l(1:floor(end/5),:)));
colormap('jet');
axis xy;
xlabel('Frequenza (Hz)');
ylabel('Frequenza (dB)');

% %Forma d'onda e Spettrogramma R
% figure('Name','Forma d^onda e Spettrogramma R','NumberTitle','off');
% subplot(2,1,1);
% plot(T_sign_r, Y_r);
% xlabel('Time');
% ylabel('Amplitude');
% subplot(2,1,2);
% imagesc(T_sp_r,F_r,mag2db(STFTmag_r));
% colormap('jet')
% axis xy;

%% SPETTRO E FASE L & R

% %Spettro e Fase L
% figure('Name','Spettro e Fase Segnale R','NumberTitle','off');
% subplot(2,1,1);
% semilogx(F_l,(diffW_l(:,25)));
% xlabel('Freq');
% hold on
% % semilogx(((sr/2)*(locPicchi_l(:,25) - 1)/length(F_l)),(valFase_l(:,25)),'x');
% xlabel('Freq');
% hold off
% subplot(2,1,2);
% semilogx(F_l,(mag2db(STFTmag_l(:,25))));
% hold off;
% xlabel('Freq');
% ylabel('dB');
% hold on
% semilogx(((sr/2)*(locPicchi_l(:,25) - 1)/length(F_l)),valPicchi_l(:,25),'x');
% xlabel('Freq');
% hold off
% 
% %Spettro e Fase R
% figure('Name','Spettro e Fase Segnale L','NumberTitle','off');
% subplot(2,1,1);
% semilogx(F_r,(diffW_r(:,25)));
% xlabel('Freq');
% hold on
% % semilogx(((sr/2)*(locPicchi_r(:,25) - 1)/length(F_r)),(valFase_r(:,25)),'x');
% xlabel('Freq');
% hold off
% subplot(2,1,2);
% semilogx(F_r,(mag2db(STFTmag_r(:,25))));
% hold off;
% xlabel('Freq');
% ylabel('dB');
% hold on
% semilogx(((sr/2)*(locPicchi_r(:,25) - 1)/length(F_r)),valPicchi_r(:,25),'x');
% xlabel('Freq');
% ylabel('dB');
% hold off
%%
subplot(2,1,1);
plot(T_sign_l,Y_l);
title("Left"); 
xlabel("Time");
ylabel("Amp");
hold on
plot(T_din_l,(picco_left),'r','linewidth',2);
hold on
plot(T_din_l,(rms_left),'black', 'linewidth',2);
hold on
plot(T_din_l,(crestF_left./max(crestF_left)));
hold off
legend('waveform','picco','rms', 'Crest factor');

subplot(2,1,2);
plot(T_sign_r,Y_r);
title("Right"); 
xlabel("Time");
ylabel("Amp");
hold on
plot(T_din_l,(picco_right),'r','linewidth',2);
hold on
plot(T_din_l,(rms_right),'black', 'linewidth',2);
hold on
plot(T_din_l,(crestF_right./max(crestF_right)));
hold off
legend('waveform','picco','rms', 'Crest factor');



%% PLOT PICCO - RMS - CREST FACTOR
 
% %Stampo i risultati : Segnale - Picco - RMS - Crest Factor
figure('Name','Signal','NumberTitle','off');
% %Stampo Left
subplot(4,2,1)                      %Left
plot(T_sign_l,Y_l, 'grey');
title("Guitar"); 
xlabel("Time");
ylabel("Amp");
% %Stampo right
subplot(4,2,2);                     %Right
% plot(T_sign_r,Y_r);
% title("Right"); 
% xlabel("Time");
% ylabel("Amp");
% %Stampo Picco
subplot(4,2,3);                     %Picco Left
plot(T_din_l,(picco_left));
title("Picco left"); 
xlabel("Time");
ylabel("Amp");
subplot(4,2,4);                     %Picco right
% plot(T_din_r,(picco_right));
% title("Picco right"); 
% xlabel("Time");
% ylabel("Amp");
% %Stampo Rms
subplot(4,2,5);                     %Rms Left
plot(T_din_l,mag2db(rms_left));
title("Rms left"); 
xlabel("Time");
ylabel("Rms");
subplot(4,2,6); 
% % Rms Right
% plot(T_din_r,mag2db((rms_right)));
% title("Rms right"); 
% xlabel("Time");
% ylabel("Rms");
% ylabel("dB");
% %Stampo Crest Factor
subplot(4,2,7);                     %Crest Factor Left
plot(T_din_l,mag2db(crestF_left));
title("Crest Factor left"); 
xlabel("Time");
ylabel("dB");
subplot(4,2,8);                     %Crest Factor Right
% plot(T_din_r,mag2db(crestF_right));
% title("Crest Factor right"); 
% xlabel("Time");
% ylabel("dB");


%% TEXT FILE EXPORT
clc
fid = fopen('Dream Theatre - Wither, Quarto Onset.txt', 'wt');

fprintf(fid, '%s\n', 'DIFFERENZE TAKE AUDIO - Dream Theatre, Quarto Onset');
fprintf(fid, '\n%s : %d\n%s : %.2f\n%s : %s\n', 'Finestra FFT', nfft, 'Overlap', round(overlap_f, 2), 'Tipo Finestra', 'Hamming');
fprintf(fid, '\n%s\n', 'DINAMICA');
fprintf(fid, '\n%s : %f dB\n%s : %f dB\n%s : %f dB\n%s : %f dB\n%s : %f dB\n%s : %f dB\n', 'Media di RMS (L - R)', (mag2db(rms_meanLeft) - mag2db(rms_meanRight)), ...
    'Deviazione std di RMS (L - R)', (mag2db(rms_stdLeft) - mag2db(rms_stdRight)), 'Media di Picco (L - R)', (mag2db(picco_meanLeft) - mag2db(picco_meanRight)), ...
    'Deviazione std di Picco (L - R)', (mag2db(picco_stdLeft) - mag2db(picco_stdRight)), 'Media di Crest Factor (L - R)', (mag2db(crestF_meanLeft) - mag2db(crestF_meanRight)), ...
    'Deviazione std di Crest Factor (L - R)', (mag2db(crestF_stdLeft) - mag2db(crestF_stdRight)));
fprintf(fid, '\n%s\n','F0 Difference & F0 Value');
fprintf(fid, '\nLeft : %f Hz\nRight : %f Hz\n', F0_lmedian, F0_rmedian);
fprintf(fid, '\n%s : %f cents\n%s : %f cents\n', 'Media di F0 (L - R)', F0mean, 'Deviazione std di F0 (L - R)', F0std);
fprintf(fid, '\n%s\n', 'Harmonic Difference');
for harm = 1 : nharm
    fprintf(fid, '\n%d %s : %f cents\n', harm, 'Armonica', PHmedianDiffHarm(harm));
end
fprintf(fid, '\n%s\n', 'Onset Delay');
for c = 1 : size(newL, 2)
    fprintf(fid, '%d %s : %f sec\n', c, 'Onset', lagD(c));
end
