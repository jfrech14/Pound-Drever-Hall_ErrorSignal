%For a pherical symmetric cavity

L=.135;  %cavity length in meters
C=2.998E8; %speed of light m/s
t=0;  %stating time is static for future use
r1=5; %radius of M1 in meters
g1=1-L/r1;
r2=10000000000000000;  %radius of M2 in meters
g2=1-L/r2;
lambda=811.754E-9; %wavelength in meters
omega0=2*pi*C/lambda;  %light angular frequency
nu0=omega0/2/pi; %light frequency
R1=0.9973;   %Mirror 1 reflectivity
R2=0.9973;   %Mirror 2 reflectivity
r=sqrt(R1*R2);
Trt=2*L/C;  %Round-trip time
linewidth=C/2/L*(1-sqrt(R1*R2))/(pi*((R1*R2)^.25)); %cavity linewidth
Finesse=(pi*(R1*R2)^(.25))/(1-sqrt(R1*R2));
FSR=C/(2*L);   %free spectral range
Q=C/linewidth/lambda; %Cavity Q-factor
Lifetime=Q/omega0; %Photon Lifetime

for m=0:1:3
    for n=0:1:3
        Modes(m+1,n+1)=C/2/L*(1+m+n)/pi*acos(sqrt(g1*g2))/1000000;
    end
end


Feom=35.15E6;  %EOM frequency
OMEGA=2*pi*Feom;  %EOM angular frequency
Pc=0.016;  %Carrier power in W
Ps=0.00072;  %Sideband power in W
Pss=0.000005;  %Secondary sideband power in W


omega=omega0-3*OMEGA:1E5:omega0+3*OMEGA;
F=r*(exp(i.*(omega-omega0)./FSR)-1)./(1-r*r*exp(i.*(omega-omega0)./FSR));
F2=F.*conj(F);
FmO=r*(exp(i.*((omega-OMEGA)-omega0)./FSR)-1)./(1-r*r*exp(i.*((omega-OMEGA)-omega0)./FSR));
FmO2=FmO.*conj(FmO);
FpO=r*(exp(i.*((omega+OMEGA)-omega0)./FSR)-1)./(1-r*r*exp(i.*((omega+OMEGA)-omega0)./FSR));
FpO2=FpO.*conj(FpO);
FmmO=r*(exp(i.*((omega-2*OMEGA)-omega0)./FSR)-1)./(1-r*r*exp(i.*((omega-2*OMEGA)-omega0)./FSR));
FmmO2=FmmO.*conj(FmmO);
FppO=r*(exp(i.*((omega+2*OMEGA)-omega0)./FSR)-1)./(1-r*r*exp(i.*((omega+2*OMEGA)-omega0)./FSR));
FppO2=FppO.*conj(FppO);

%Only includes 1st order sidebands (Pinc for 2nd order is not calculated)
Pincident=Pc*F2+Ps*(FpO2+FmO2)+2*sqrt(Pc*Ps)*real(F.*conj(FpO)-...
    conj(F).*FmO)*cos(OMEGA*t)+2*sqrt(Pc*Ps)*imag(F.*conj(FpO)-...
    conj(F).*FmO)*sin(OMEGA*t);
error1=-2*sqrt(Pc*Ps)*imag(F.*conj(FpO)-conj(F).*FmO);

%Error signal with 1st & 2nd order sidebands for higher EOM powers
error2=-2*sqrt(Pc*Ps)*imag(F.*conj(FpO)-conj(F).*FmO)...
    -2*sqrt(Ps*Pss)*imag(FpO.*conj(FppO)-conj(FpO).*FmmO)-2*sqrt(Ps*Pss)*imag(FpO.*conj(FmmO)-conj(FpO).*FppO) -2*sqrt(Ps*Pss)*imag(FmO.*conj(FppO)-conj(FmO).*FmmO)-2*sqrt(Ps*Pss)*imag(FmO.*conj(FmmO)-conj(FmO).*FppO);

angfrequency=-3*OMEGA/(10^6):0.1:3*OMEGA/(10^6);
frequency=angfrequency/2/pi;
plot(frequency,error1)