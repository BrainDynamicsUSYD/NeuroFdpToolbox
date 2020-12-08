%% This function plots the power spectrum of the signals using FFT in the 
% loglog plot and also rescale x-axis to be equally distributed
%
% author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% dataFileName = 'my144_101' ;
dataFileName = 'my147_02' ;
% dataFileName = 'ma027_032' ;
% dataFileName = 'ma025_03' ;
cd ..
cd ..
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigIn,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
clearvars LFPs

%% Average Power Spectrum of the signals
% choose 1 spatial point
xPoint = 5 ;
yPoint = 5 ;
% You can compute the Fourier transform of an image in Matlab using fft2 as follows:
nfft = 10000 ;
%dataNew = dataFull(2:8,2:8,:) ;
%dataSelect = dataNew(:,:,round(100*fsTemporal): round(180*fsTemporal) ) ;
dataSelect = squeeze(sigIn(xPoint,yPoint,ceil(1*fsTemporal):ceil(300*fsTemporal)) );

vf = fft((dataSelect(:)),nfft) ;
%vf2 = fft(squeeze(filteredLFPs(xPoint,yPoint,:)),nfft) ;
% The function fftshift is needed to put the DC component (frequency = 0) in the center.
%The power spectrum is just the square of the modulus of the Fourier transform,  which is obtained as follows:
vpf = abs(vf).^2 ;
%vpf2 = abs(vf2).^2 ;
% To display the power spectrum you will need to define some frequency coordinates.
% If the image is of size N, then the frequencies run from -N/2 to N/2-1 (assuming N is even):
% [~,~,N] = size(LFPs) ;
freqIndex = 1 : (nfft/2+1) ;
freqAxis = freqIndex/nfft*fsTemporal ;
nPts = 400 ;
% fnew=fsTemporal/2.*logspace(log10(fsTemporal/length(vpf)),0,nPts);
fnew=fsTemporal/2.*logspace(log10(fsTemporal/length(vpf)/fsTemporal),0,nPts);
Ynew= interp1(freqAxis,vpf(1:nfft/2+1),fnew);

% Then display the log power spectrum using imagesc:
%close all
figure; set(gcf,'Visible','On')
loglog(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
% loglog(fnew(100:370),Ynew(100:370), 'LineWidth', 1.6);
xlim([1 75])
title('Power spectrum of the electrode [5,5] over 300s')
ax.TitleFontSizeMultiplier = 4;
ylabel('Power')
xlabel('Frequency(Hz)')
% figure; set(gcf,'Visible','On')
% plot(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
% loglog(freqAxis,vpf2(freqIndex));
% xlim([1 500])


%%
freqIndex = 1 : (nfft/2+1) ;
freqAxis = freqIndex/nfft*fsTemporal ;
nPts = 400 ;
% fnew=fsTemporal/2.*logspace(log10(fsTemporal/length(vpf)),0,nPts);
fnew=fsTemporal/2.*logspace(log10(fsTemporal/length(vpf)/fsTemporal),0,nPts);
Ynew= interp1(freqAxis,vpf(1:nfft/2+1),fnew);

figure; set(gcf,'Visible','On')
% loglog(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
loglog(fnew(100:370),Ynew(100:370), 'LineWidth', 1.6);
loglog(fnew,Ynew, 'LineWidth', 1.6);
xlim([0.1 500])
title('Power spectrum of the electrode [5,5] over 300s')
ax.TitleFontSizeMultiplier = 4;
ylabel('Power')
xlabel('Frequency(Hz)')






