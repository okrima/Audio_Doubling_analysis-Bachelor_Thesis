
% Alcuni parametri
spectralThres = 0.5;
% hop = 512;
% w = hann(1024);
% 
% % Leggo un file (mono, è meglio)
% [x,sr] = audioread('64.wav');
% x = mean(x,2);
% 
% % STFT
% STFT = getFD(x, sr, hop, w);

% scalo tra 0 e 1, con -96dB come soglia inferiore
M = amp2db(abs(STFT_l), -96) - 96;
M = M./max(M(:));

% Percussive mask
D = (20 * log10(M(:,1:end-1)./M(:,2:end))) > spectralThres;

% Percussive feature
feature = sum(D, 1);

figure('Name', 'Presti');
plot(T_sp_l(1:end-1), feature);
hold on
plot(T_sign_l, Y_l*1000);
