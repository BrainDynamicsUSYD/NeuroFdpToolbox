% log bin fft with interpolation
y = Gamma_one;
npts = 20 ;
NFFT = length(Gamma_one) ;
Y = abs(fft(Gamma_one)) ;
Y = Y(1:NFFT/2+1) ;
f = linspace(0,fsTemporal/2,NFFT/2+1) ;

fnew=fsTemporal/2.*logspace(0,log10(fsTemporal/length(y)),npts);
Ynew= interp1(f,Y(1:NFFT/2+1),fnew);
figure;
loglog(fnew,Ynew,'.','markerSize',24)