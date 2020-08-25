function LSD_2(Clean,Noisy,fs,RS,Range)

%LSD=LOGSPECTRALDISTANCE(CLEAN,NOISY,FS,RS,RANGE)
% Calculates the average log-spectral distance between CLEAN and NOISY
% signals. Frames are 25ms long with 60 percent (15ms) overlap, hamming
% windowed. RS is the remove silence option (default: 0) if 1, the program
% uses a Voice Activity Detector to find non-speech frames and eliminates
% them from calculation of LSD. FS is the sampling frequency (Default
% 10KHz). RANGE is the frequency range used in the calculation of LSD
% (default: [0 FS/2]) it is a two-element vector the first element should
% be smaller than the second one and non of them should be greater than
% FS/2.

if nargin<3
    fs=10000;
end
if nargin<4
    RS=0;
end
if nargin<5
    Range=[0,fs/2];
end
if RS==1
    [Clean, Noisy]=RemoveSilence(Clean,Noisy,fs);
end

LPO=round(24*fs/16000);
Nf=256; %Frequency Responce Length

Len=min(length(Clean),length(Noisy));
Clean=Clean(1:Len);
Noisy=Noisy(1:Len);
Clean=Clean./sqrt(sum(Clean.^2));
Noisy=Noisy./sqrt(sum(Noisy.^2));

W=round(.025*fs);
SP=.4;

CL=segment(Clean,W,SP);
NO=segment(Noisy,W,SP);

[Ac, Gc]=lpc(CL,LPO);
[An, Gn]=lpc(NO,LPO);

RW=2*Range*pi*(1-1/(Nf+1))/fs;
w=linspace(RW(1),RW(2),Nf);
WFZ=zeros(LPO+1,Nf);
IMAGUNIT=sqrt(-1);
for k=0:LPO
    WFZ(k+1,:)=exp(IMAGUNIT*k*w);
end

nframes=size(Ac,1);
LSD=zeros(1,nframes);
for i=1:nframes
    Sc=sqrt(Gc(i))./(abs(Ac(i,:)*WFZ));
    Sn=sqrt(Gn(i))./(abs(An(i,:)*WFZ));
    LSD(i)=sqrt(mean((log(Sc)-log(Sn)).^2));
end

LSD=mean(LSD);
