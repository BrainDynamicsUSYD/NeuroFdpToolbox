function powerSpectrumFFT(sigIn,Fs)
% This function finds the power spectrum by FFT
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% choose 1 spatial point
xPoint = 1 ;
yPoint = 1 ;
% You can compute the Fourier transform of an image in Matlab using fft2 as follows:
nfft = 4000 ;
vf = fft(squeeze(sigIn(xPoint,yPoint,:)),nfft) ;
%vf2 = fft(squeeze(filteredLFPs(xPoint,yPoint,:)),nfft) ;
% The function fftshift is needed to put the DC component (frequency = 0) in the center.  
%The power spectrum is just the square of the modulus of the Fourier transform,  which is obtained as follows:
vpf = abs(vf).^2 ;
%vpf2 = abs(vf2).^2 ;
% To display the power spectrum you will need to define some frequency coordinates.  
% If the image is of size N, then the frequencies run from -N/2 to N/2-1 (assuming N is even):
% [~,~,N] = size(LFPs) ;
freqIndex = 1 : ceil(nfft/2) ;
freqAxis = freqIndex/nfft*Fs ;
% Then display the log power spectrum using imagesc:
% close all
figure; set(gcf,'Visible','On')
% plot(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
loglog(freqAxis,vpf(freqIndex)); 
xlim([1 400])
xlabel('frequency (Hz)')
ylabel('power')
title('power spectrum')
%figure; set(gcf,'Visible','On')
% plot(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
%loglog(freqAxis,vpf2(freqIndex)); 
%xlim([1 500])