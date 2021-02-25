clear;
tic
c=3*10^8;
radius=0.305/2;
% radius=0.127/2;
freq = 2.40E9:0.3125E6:2.5E9;
% freq = 59.4E9:2E7:60.4E9;
% freq = 2.40E9:0.2E7:2.5E9;

% besDerZerMatD = load('besDerZerMat5k.mat');
% besDerZerMat = (besDerZerMatD.besDerZerMat)';
% besZerMatD = load('besZerMat5k.mat');
% besZerMat = (besZerMatD.besZerMat)';

load('besDerZerMat.mat');
load('besZerMat.mat');

Zo=50;
eta = 377;
l = 0.035;
% l = 0.005;

sigma = 1E6;
mu = 4*pi*10^-7;
R = sqrt((2*pi*freq*mu)./(2*sigma));
k = 2*pi*freq./c;
% WGlen = 4;
WGlen = 4.57;
% WGlen = 1;
channel = [];
att = [];

nTE=[];
mTE=[];
fcTE=[];
coWnTE=[];
nTM=[];
mTM=[];
fcTM=[];
coWnTM=[];
for fi=1:length(freq)
    for m=1:1001
        for n=1:1000
            fcTEtemp=(c/(2*pi*radius))*besDerZerMat(m,n);
            fcTMtemp=(c/(2*pi*radius))*besZerMat(m,n);
            if fcTEtemp <= freq(fi) % && fcTEtemp >= 57E9 % 1.5E9
                nTE = [nTE; n];
                mTE = [mTE; m-1];
                fcTE = [fcTE; fcTEtemp];
                coWnTE = [coWnTE; besDerZerMat(m,n)/radius];
            end
            if fcTMtemp <= freq(fi) % && fcTMtemp >= 57E9 % 1.5E9
                nTM = [nTM; n];
                mTM = [mTM; m-1];
                fcTM = [fcTM; fcTMtemp];
                coWnTM = [coWnTM; (besZerMat(m,n))/radius];
            end
        end
    end
end
[radresTE, radresTM, gammaTE, gammaTM] = radResCyl_multitone(mTE,nTE,mTM,nTM,radius,freq,fcTE,fcTM,c,k,R,eta,l,coWnTE,coWnTM);

radreacTE = imag(hilbert(radresTE));
radreacTM = imag(hilbert(radresTM));

for fi=1:length(freq)
    antimpTE(fi) = sum(radresTE(fi,:)) + 1i*sum(radreacTE(fi,:));
    antimpTM(fi) = sum(radresTM(fi,:)) + 1i*sum(radreacTM(fi,:));
    for n = 1:length(nTE)
        TEmodeimp(fi,n)=radresTE(fi,n)+1i*radreacTE(fi,n);
    end
    for n = 1:length(nTM)
        TMmodeimp(fi,n)=radresTM(fi,n)+1i*radreacTM(fi,n);
    end    
    TsTE = diag(exp(-1*gammaTE(fi,:)*WGlen));
    chTEmode = TEmodeimp(fi,:)*TsTE;
    TsTM = diag(exp(-1*gammaTM(fi,:)*WGlen));
    chTMmode = TMmodeimp(fi,:)*TsTM;
    channel(fi) = ((2*Zo)/(abs(antimpTM(fi) + Zo + antimpTE(fi))^2))*...
        (sum(chTEmode)+sum(chTMmode));
%     if isnan(channel(fi)) == 1
%         channel(fi) = 0;
%         att(fi) = 0;
%     else
%         att(fi) = 10*log10(abs(channel(fi)));
%     end
    att(fi) = 10*log10((abs(channel(fi)))^2);
%     att(fi) = 10*log10(abs(channel(fi)));
end

figure
plot(freq,att);
toc
