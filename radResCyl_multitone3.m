function [radresTE, radresTM, gammaTE, gammaTM] = radResCyl_multitone3(m_TE,n_TE,m_TM,n_TM,radius,freq,fcTE,fcTM,k,Rs,eta,l,coWnTE,coWnTM)
nmTE(:,1) = m_TE;
nmTE(:,4) = coWnTE/radius;
nmTE(:,2) = n_TE;
nmTM(:,1) = m_TM;
nmTM(:,2) = n_TM;
nmTM(:,4) = coWnTM/radius;
nmTE(:,3) = fcTE;
nmTM(:,3) = fcTM;
for fi = 1:length(freq)
    for n = 1:length(nmTE(:,1))
        pnmSquare = (nmTE(n,4))^2;
        nSquare = (nmTE(n,1))^2;
        kSquare = (k(fi))^2;
        besselTemp = (besselj(nmTE(n,1),nmTE(n,4))) * ...
            (besselj(nmTE(n,1),nmTE(n,4)));
        sinSqu = (sin(l*k(fi)*pi/180))*(sin(l*k(fi)*pi/180));
        fcByFreqTESq = (nmTE(n,3)/freq(fi))^2;
        lBiRadius = l/radius;
        TEMat(n,1,fi) = (sqrt(kSquare - pnmSquare*radius*radius));
        alphaTempTE = (Rs(fi)*k(fi))/(eta*radius*TEMat(n,1,fi));
        TEMat(n,2,fi) = 8.85*alphaTempTE*(fcByFreqTESq + (nSquare/(pnmSquare - nSquare)));
        gammaTE(fi,n) = TEMat(n,2,fi) + 1i*TEMat(n,1,fi);
        constRadTE = eta*k(fi)*nSquare/(TEMat(n,1,fi)*pi*sinSqu*...
            (pnmSquare - nSquare)*besselTemp);
        radTEfunc=@(zeta) besselj(nmTE(n,1),zeta*nmTE(n,4)).*...
            sin(k(fi)*radius*(lBiRadius-1+zeta)*pi/180)./zeta;
        inteTempTE = integral(radTEfunc,1-lBiRadius,1);
        TEMat(n,3,fi) = constRadTE*inteTempTE*inteTempTE;
        radresTE(fi,n) = TEMat(n,3,fi);
    end
    for n = 1:length(nmTM(:,1))
        pnmSquare = (nmTM(n,4))*(nmTM(n,4));
        sinSqu = (sin(l*k(fi)*pi/180))*(sin(l*k(fi)*pi/180));
        lBiRadius = l/radius;
        TMMat(n,1,fi) = (sqrt(k(fi)^2 - pnmSquare*radius*radius));
        alphaTempTM = (Rs(fi)*k(fi))/(eta*radius*TMMat(n,1,fi));
        TMMat(n,2,fi) = 8.85*alphaTempTM;
        gammaTM(fi,n) = TMMat(n,2,fi) + 1i*TMMat(n,1,fi);
        besselDeri = besselj(nmTM(n,1)-1,nmTM(n,4)) - ...
        (nmTM(n,1)/nmTM(n,4))*besselj(nmTM(n,1),nmTM(n,4));
        besselTemp = besselDeri * besselDeri;
        constRadTM = eta*TMMat(n,1,fi)/(k(fi)*pi*sinSqu*besselTemp);
        radTMfunc=@(zeta) (besselj(nmTM(n,1)-1,zeta*nmTM(n,4)) - (nmTM(n,1)./...
            (zeta*nmTM(n,4))).*besselj(nmTM(n,1),zeta*nmTM(n,4))).*...
            sin(k(fi)*radius*(lBiRadius-1+zeta)*pi/180);
        inteTempTM = integral(radTMfunc,1-lBiRadius,1);
        TMMat(n,3,fi)= constRadTM*inteTempTM*inteTempTM;
        radresTM(fi,n) = TMMat(n,3,fi);
    end
end
end

