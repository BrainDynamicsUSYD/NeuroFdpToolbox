function job_spectrum(arrayID)

cd ..
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%
% arrayID = 1 ;

switch arrayID
    case 1
        dataFileName = 'ma027_01' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')

    case 2
        dataFileName = 'ma027_02' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')

    case 3
        dataFileName = 'ma027_03' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')

    case 4
        dataFileName = 'ma027_04' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
        
    case 5
        dataFileName = 'ma027_05' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
        
    case 6
        dataFileName = 'my147_01' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
        
    case 7
        dataFileName = 'my147_02' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')

    case 8
        dataFileName = 'my147_03' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
        
    case 9
        dataFileName = 'my147_04' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
        
    case 10
        dataFileName = 'my147_05' ; 
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')
    
    case 11    
        dataFileName = 'ma027_032' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'Fs')

    case 12    
        dataFileName = 'my147_53' ;
        load([pwd,'/Data/UtahArrayData/',dataFileName],'sampleFs')
        Fs = sampleFs.LFP ;

end


%% load data
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs')

%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
flagBandstop = 0 ;
[sigOriRaw,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
sigOriRS = reshape(sigOri, 100,[]) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% power spectrum (average)
% nfft = size(sigOriRS,2) ;
% nfft = 10000 ;
% vpf = zeros(1,nfft) ;
% 
% for iChannel = setdiff(1:100, badChannels)
%     vf = fft(squeeze(sigOriRS(iChannel,:)),nfft) ;
%     % The function fftshift is needed to put the DC component (frequency = 0) in the center.
%     %The power spectrum is just the square of the modulus of the Fourier transform,  which is obtained as follows:
%     vpf = vpf + abs(vf).^2 ; 
% end
% vpf = vpf / iChannel;
% 
% % To display the power spectrum you will need to define some frequency coordinates.
% % If the image is of size N, then the frequencies run from -N/2 to N/2-1 (assuming N is even):
% % [~,~,N] = size(LFPs) ;
% freqIndex = 1 : ceil(nfft/2) ;
% freqAxis = freqIndex/nfft*Fs ;
% % Then display the log power spectrum using imagesc:
% % close all
% figure; set(gcf,'Visible','On')
% % plot(freqAxis,vpf(freqIndex))
% %subplot(1,2,1)
% loglog(freqAxis,vpf(freqIndex));
% xlim([0.5 400])
% xlabel('frequency (Hz)')
% ylabel('power')
% title(['Average spectrum of ', dataFileName(1:5),'\_',dataFileName(8)...
%     ' over ',num2str(fix(size(sigOriRS,2)/Fs)),'s'])
% %figure; set(gcf,'Visible','On')
% % plot(freqAxis,vpf(freqIndex))
% %subplot(1,2,1)
% %loglog(freqAxis,vpf2(freqIndex));
% %xlim([1 500])
% t = datetime('now') ;
% dateStr = datestr(t,'mmmmdd_HH:MM') ;
% savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'PowerSpectrum',dateStr,'.fig'])
% 
% %% power spectrum (10s)
% timeStep = 2 ;
% timeEpoch = 10000 ;
% timeStart = timeEpoch*(timeStep-1)  ;
% timeEnd = timeEpoch*timeStep   ;
% 
% sigOriRS = reshape(sigOri(:,:,timeStart:timeEnd), 100,[]) ;
% nfft = size(sigOriRS,2) ;
% % nfft = 10000 ;
% vpf = zeros(1,nfft) ;
% 
% for iChannel = 55 % setdiff(1:100, badChannels)
%     vf = fft(squeeze(sigOriRS(iChannel,:)),nfft) ;
%     % The function fftshift is needed to put the DC component (frequency = 0) in the center.
%     %The power spectrum is just the square of the modulus of the Fourier transform,  which is obtained as follows:
%     vpf = vpf + abs(vf).^2 ; 
% end
% % vpf = vpf / iChannel;
% 
% % To display the power spectrum you will need to define some frequency coordinates.
% % If the image is of size N, then the frequencies run from -N/2 to N/2-1 (assuming N is even):
% % [~,~,N] = size(LFPs) ;
% freqIndex = 1 : ceil(nfft/2) ;
% freqAxis = freqIndex/nfft*Fs ;
% 
% %figure; set(gcf,'Visible','On')
% % plot(freqAxis,vpf(freqIndex))
% %subplot(1,2,1)
% %loglog(freqAxis,vpf2(freqIndex));
% %xlim([1 500])
% t = datetime('now') ;
% dateStr = datestr(t,'mmmmdd_HH:MM') ;
% savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'PowerSpectrum10s',dateStr,'fig'])
% figure; set(gcf,'Visible','On')
% % plot(freqAxis,vpf(freqIndex))
% %subplot(1,2,1)
% loglog(freqAxis,vpf(freqIndex));
% xlim([0.5 400])
% xlabel('frequency (Hz)')
% ylabel('power')
% title(['10s spectrum of ', dataFileName(1:5),'\_',dataFileName(8)...
%     ' over ',num2str(fix(size(sigOriRS,2)/Fs)),'s'])
% %figure; set(gcf,'Visible','On')
% % plot(freqAxis,vpf(freqIndex))
% %subplot(1,2,1)
% %loglog(freqAxis,vpf2(freqIndex));
% %xlim([1 500])
% 
% savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'PowerSpectrum',dateStr,'.fig'])
% 
% %% wavelet spectrogram (10s)
% timeStep = 2 ;
%     addpath('ToolOthers/uimage')
%     timeEpoch = 10000 ;
%     timeStart = timeEpoch*(timeStep-1)  ;
%     timeEnd = timeEpoch*timeStep   ;
%     %timeStart = 0 ;
%     %timeEnd = 20 ;
%     
%     freqRange = 0.5:100 ;
%     fc = centfrq('cmor1.5-1') ;
%     scalerange = fc./(freqRange/fsTemporal) ;
%     scales = scalerange(end):0.5:scalerange(1) ;
%     pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
%     
%     tempData =  squeeze (sigOri(5,5,timeStart:timeEnd)) ; 
%     
%     wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
%     tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
%     figure
%     % wtNorm = zscore(abs(wt),[],2) ;
%     
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,:)) )))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'NormWaveletSpectrogram10s',dateStr,'.fig'])
% 
%     figure
%     % uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'WaveletSpectrogram10s',dateStr,'.fig'])
% 
% 
% %% wavelet spectrogram (2s)
% timeStep = 2 ;
%     addpath('ToolOthers/uimage')
%     timeEpoch = 2000 ;
%     timeStart = timeEpoch*(timeStep-1)  ;
%     timeEnd = timeEpoch*timeStep   ;
%     %timeStart = 0 ;
%     %timeEnd = 20 ;
%     
%     freqRange = 0.5:100 ;
%     fc = centfrq('cmor1.5-1') ;
%     scalerange = fc./(freqRange/fsTemporal) ;
%     scales = scalerange(end):0.5:scalerange(1) ;
%     pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
%     
%     tempData =  squeeze (sigOri(5,5,timeStart:timeEnd)) ; 
%     
%     wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
%     tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
%     figure
%     % wtNorm = zscore(abs(wt),[],2) ;
%     
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,:)) )))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'NormWaveletSpectrogram2s',dateStr,'.fig'])
% 
%     figure
%     % uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'WaveletSpectrogram2s',dateStr,'.fig'])
% 
% %% wavelet spectrogram (10s)
% timeStep = 2 ;
%     addpath('ToolOthers/uimage')
%     timeEpoch = 2000 ;
%     timeStart = timeEpoch*(timeStep-1)  ;
%     timeEnd = timeEpoch*timeStep   ;
%     %timeStart = 0 ;
%     %timeEnd = 20 ;
%     
%     freqRange = 20:100 ;
%     fc = centfrq('cmor1.5-1') ;
%     scalerange = fc./(freqRange/fsTemporal) ;
%     scales = scalerange(end):0.5:scalerange(1) ;
%     pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
%     
%     tempData =  squeeze (sigOri(5,5,timeStart:timeEnd)) ; 
%     
%     wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
%     tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
%     figure
%     % wtNorm = zscore(abs(wt),[],2) ;
%     
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,:)) )))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'NormWaveletSpectrogramHighFreq',dateStr,'.fig'])
% 
%     figure
%     % uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'WaveletSpectrogramHighFreq',dateStr,'.fig'])
% %% wavelet spectrogram (10s)
% timeStep = 2 ;
%     addpath('ToolOthers/uimage')
%     timeEpoch = 2000 ;
%     timeStart = timeEpoch*(timeStep-1)  ;
%     timeEnd = timeEpoch*timeStep   ;
%     %timeStart = 0 ;
%     %timeEnd = 20 ;
%     
%     freqRange = 0.5:20 ;
%     fc = centfrq('cmor1.5-1') ;
%     scalerange = fc./(freqRange/fsTemporal) ;
%     scales = scalerange(end):0.5:scalerange(1) ;
%     pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
%     
%     tempData =  squeeze (sigOri(5,5,timeStart:timeEnd)) ; 
%     
%     wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
%     tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
%     figure
%     % wtNorm = zscore(abs(wt),[],2) ;
%     
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,:)) )))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'NormWaveletSpectrogramLowFreq',dateStr,'.fig'])
% 
%     figure
%     % uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
%     uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
%     savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
%     'WaveletSpectrogramLowFreq',dateStr,'.fig'])

%% multitaper analysis (10s)
PxxMean = zeros(size(sigOriRS,3),1) ;
for iChannel = setdiff(1:10,badChannels)
    [Pxx,F] = pmtm(sigOriRS(iChannel,:),[],10000,Fs) ;
    PxxMean = Pxx + PxxMean ;
end
PxxMean = PxxMean/iChannel ;
plot(F,PxxMean)
xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0.5 400])
title(['Average spectrum of ', dataFileName(1:5),'\_',dataFileName(8)...
    ' over ',num2str(fix(size(sigOriRS,2)/Fs)),'s'])

t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;
savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
    'MultitaperSpectrum',dateStr,'.fig'])

%% multitaper analysis (10s)
timeStep = 2 ;
timeEpoch = 10000 ;
timeStart = timeEpoch*(timeStep-1)  ;
timeEnd = timeEpoch*timeStep   ;

[Pxx,F] = pmtm(sigOriRS(55,timeStart:timeEnd),[],[],Fs) ;
loglog(F,Pxx)
xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0.5 400])
title(['10s spectrum of ', dataFileName(1:5),'\_',dataFileName(8)...
    ' over ',num2str(fix(size(sigOriRS,2)/Fs)),'s'])

savefig([pwd,'/Results/Project2/Results/powerspectrum/',dataFileName,...
    'MultitaperSpectrum10s',dateStr,'.fig'])
