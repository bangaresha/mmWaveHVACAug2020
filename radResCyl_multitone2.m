function [radresTE, radresTM, gammaTE, gammaTM,powerf] = radResCyl_multitone2(nm_TE,nm_TM,radius,freq,c,k,Rs,eta,l)
gammaTE = zeros(length(freq), length(nm_TE(:,1)));
gammaTM = zeros(length(freq), length(nm_TM(:,1)));
radresTE = zeros(length(freq), length(nm_TE(:,1)));
radresTM = zeros(length(freq), length(nm_TM(:,1)));
pTE_f = [];
fcTE_f = [];
mnTE_f = [];
pTM_f = [];
fcTM_f = [];
mnTM_f = [];

for n = 1:length(nm_TE)
    pnmSquare = (nm_TE(n,4))*(nm_TE(n,4));
    betaTE = (sqrt(k^2 - pnmSquare*radius*radius));
    alphaTempTE = (Rs*k)/(eta*radius*betaTE);
    nSquare = nm_TE(n,1)^2;
    fcByFreqTESq = (nm_TE(n,3)/freq)^2;
    alphaTE = 8.68*alphaTempTE*(fcByFreqTESq + (nSquare/(pnmSquare - nSquare)));
    gammaTE(n)=alphaTE+1i*betaTE;
    modeVelocityTE = sqrt(1 - fcByFreqTESq);
    besselTemp = besselj(nm_TE(n,1),nm_TE(n,4)) * besselj(nm_TE(n,1),nm_TE(n,4));
    sinSqu = (sin(l*k*pi/180))*(sin(l*k*pi/180));
    constRadTE = eta*k*nSquare/(betaTE*pi*sinSqu*(pnmSquare - nSquare)*besselTemp);
    lBiRadius = l/radius;
    radTEfunc=@(zeta) besselj(nm_TE(n,1),zeta*nm_TE(n,4)).*...
        sin(k*radius*(lBiRadius-1+zeta)*pi/180)./zeta;
    inteTempTE = integral(radTEfunc,1-lBiRadius,1);
    radresTE(n)= 2*constRadTE*inteTempTE*inteTempTE;
end
for n = 1:length(nm_TM)
    pnmSquare = (nm_TM(n,4))*(nm_TM(n,4));
    betaTM = (sqrt(k^2 - pnmSquare*radius*radius));
    alphaTempTM = (Rs*k)/(eta*radius*betaTM);
    nSquare = nm_TM(n,1)^2;
    fcByFreqTMSq = (nm_TM(n,3)/freq)^2;
    alphaTM = 8.68*alphaTempTM;
    gammaTM(n)=alphaTM+1i*betaTM;
    modeVelocityTM = sqrt(1 - fcByFreqTMSq);
    besselDeri = besselj(nm_TM(n,1)-1,nm_TM(n,4)) - ...
        (nm_TM(n,1)/nm_TM(n,4))*besselj(nm_TM(n,1),nm_TM(n,4));
    besselTemp = besselDeri * besselDeri;
    sinSqu = (sin(l*k*pi/180))*(sin(l*k*pi/180));
    constRadTM = eta*betaTM/(k*pi*sinSqu*besselTemp);
    lBiRadius = l/radius;
    radTMfunc=@(zeta) (besselj(nm_TM(n,1)-1,zeta*nm_TM(n,4)) - ...
        (nm_TM(n,1)./(zeta*nm_TM(n,4))).*besselj(nm_TM(n,1),zeta*nm_TM(n,4))).*...
        sin(k*radius*(lBiRadius-1+zeta)*pi/180);
    inteTempTM = integral(radTMfunc,1-lBiRadius,1);
    radresTM(n)= 2*constRadTM*inteTempTM*inteTempTM;
end

powerf = [];
sumRadResTE = sum(radresTE);
for i = 1:length(radresTE)
    pTE_f = [pTE_f; 100*radresTE(i)/sumRadResTE];
    mnTE_f = [mnTE_f; strcat('TE',string(nm_TE(i,1)),',' ,string(nm_TE(i,2)))];
    fcTE_f = [fcTE_f; nm_TE(i,3)];
    if (100*radresTE(i)/sumRadResTE) >= 5
        powerf = [powerf;strcat('TE',string(nm_TE(i,1)),',' ,string(nm_TE(i,2)))...
            100*radresTE(i)/sumRadResTE radresTE(i)];
    end
end

sumRadResTM = sum(radresTM);
for i = 1:length(radresTM)
    pTM_f = [pTM_f; 100*radresTM(i)/sumRadResTM];
    fcTM_f = [fcTM_f; nm_TM(i,3)];
    mnTM_f = [mnTM_f; strcat('TM',string(nm_TM(i,1)), ',', string(nm_TM(i,2)))];
    if (100*radresTM(i)/sumRadResTM) >= 5
        powerf = [powerf;strcat('TM',string(nm_TM(i,1)),',' ,string(nm_TM(i,2)))...
            100*radresTM(i)/sumRadResTM radresTM(i)];
    end
end

stem([fcTE_f', fcTM_f'],[pTE_f', pTM_f']');
title('Mode Power(%Total Power) Versus Mode Cutoff Frequency');

end

