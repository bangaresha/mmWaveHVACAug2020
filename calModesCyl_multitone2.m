besDerZerMatD = load('besDerZerMat5k.mat');
besDerZerMat = besDerZerMatD.besDerZerMat';
besZerMatD = load('besZerMat5k.mat');
besZerMat = besZerMatD.besZerMat';
c=3*10^8;
radius=0.305/2; 
% radius=0.127/2;
% freq = 2.40E9:0.3125E6:2.5E9;         
% freq = 57.24E9:2E6:59.4E9;
frequency = 2.5E9;
% frequency = 60E9;
eta = 377;
l = 0.035;
% l=0.005;
sigma = 10^6;
mu = 4*pi*1E-7;
R = ((2*pi*frequency*mu)/(2*sigma))^0.5;
k = 2*pi*frequency/c;
nm_TE = [];
nm_TM = [];
for m=1:5000
    for n=1:4999
        fc_TE_temp=(c/(2*pi*radius))*besDerZerMat(n+1,m);
        fc_TM_temp=(c/(2*pi*radius))*besZerMat(n,m);
        if fc_TE_temp <= frequency %&& fc_TE_temp >= 1.5E9 % 57E9
            nm_TE = [nm_TE;n m fc_TE_temp besDerZerMat(n+1,m)];
        end
        if fc_TM_temp <= frequency % && fc_TM_temp >= 1.5E9 % 57E9
            nm_TM = [nm_TM;n-1 m fc_TM_temp besZerMat(n,m)];
        end
    end
end

[radresTES, radresTMS, gammaTES, gammaTMS, powerf] = radResCyl_multitone2(nm_TE,...
    nm_TM,radius,frequency,c,k,R,eta,l);
