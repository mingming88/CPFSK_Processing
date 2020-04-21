clear all; close all;

%% Parameter settings
BaudRate = 16e3;
Tb = 1/BaudRate; % period
OSR = 10;
fs = OSR * 1/Tb; 

% G(t)
gt_L = 4;
tvector = 0:1/fs:gt_L*Tb;
X = 2/Tb*tvector-4; %normalised
k = 1;
gt = k*sinc(X)./(1-X.^2);  
gt(find((1-X.^2)==0)) = 0.5; 
% figure(1); subplot(211); plot(tvector/Tb, gt);

Int_gt = cumsum(gt)*1/fs;
gt = gt./Int_gt(end)*0.5; %Scaling by half
% figure(1); subplot(212); plot(tvector/Tb, gt);

%% CPFSK Modulation with AWGN noise

% sig 1
Nsym = 48*20;
bitSeq = sign(-1+(2).*rand(1,Nsym));
phaseVec = cumsum(conv(upsample(bitSeq,OSR),gt))*1/fs;
foffset_norm = 0;
h = 0.5; 
sig1 = exp(1i*2*pi*h*phaseVec).*exp(1i*2*pi*foffset_norm*(0:length(phaseVec)-1));
sigpwr = mean(sig1.*conj(sig1));
sig1_unitpwr = sig1./sigpwr;
awgnnoise = 1/sqrt(2)*(randn(1,length(sig1))+1i*randn(1,length(sig1)));
desiredSNR = 10;
addednoise = 10^(-desiredSNR/20)*awgnnoise; %10^(-desiredSNR/20)*sqrt(OSR)*awgnnoise;

% sig 2
sig2_unitpwr = sig1_unitpwr;
% Nsym = 48*20;
% bitSeq = sign(-1+(2).*rand(1,Nsym));
% phaseVec = cumsum(conv(upsample(bitSeq,OSR),gt))*1/fs;
% foffset_norm = 0.1;
% h = 0.5; 
% sig2 = exp(1i*2*pi*h*phaseVec).*exp(1i*2*pi*foffset_norm*(0:length(phaseVec)-1));
% sigpwr = mean(sig2.*conj(sig2));
% sig2_unitpwr = sig2./sigpwr;

samplesshift = 12;
sig2_shift = [zeros(1,samplesshift) sig2_unitpwr(1:end-samplesshift)];

combinedsig = sig1_unitpwr+sig2_shift+addednoise;
combinedsig = [zeros(1,10) combinedsig zeros(1,samplesshift*10)];
figure; plot(fftshift(makeFreq(length(combinedsig),fs)),20*log10(fftshift(abs(fft(combinedsig)))));

sig_cutout = sig1_unitpwr;
shifts = 1:samplesshift*2; % the preamble must fall in the first 9 seconds
pri_cumu = cumsum(combinedsig.*conj(combinedsig));
[pp, f_inds] = dot_xcorr_singleChannel_2018_24(conj(sig_cutout), combinedsig, pri_cumu, norm(sig_cutout)^2, shifts);
figure; plot(shifts/fs,pp);




% CM20(sig,fs,1);

