function preAnalysisP2(wvcfs,sigOri,sigSmooth,fsTemporal)

% figure;
% colorArray = [{'g'},{'r'},{'c'},{'m'},{'k'},{'k'}] ;
% count = 1 ;
% timeRange = fix(70*fsTemporal)+200:10:fix(74*fsTemporal)+200 ;
% plot(1:length(timeRange),zscore(squeeze(sigSmooth(5,5,timeRange)))-3*3,'b','LineWidth',1)
% hold on
% % plot(1:length(timeRange),zeros(length(tempPlot))-3*3,'k' ,'LineWidth',1)
% % hold on
% 
% for iBand = 5:2:13
%     temp = (wvcfs{iBand}(5,5,:)) ;
%     prcPara = 95 ;
%     prcBound(iBand) = prctile(abs(temp),prcPara) ;
%     [tempPlot,mu,sigma] = zscore(real(temp(timeRange))) ;
%     prcBound(iBand) = (prcBound(iBand))/sigma;
%     plot(1:length(tempPlot),(squeeze(real(tempPlot))) -iBand*3,[colorArray{count},'-'],'LineWidth',1)
%     hold on
%     plot(1:length(tempPlot),(abs(squeeze(temp(timeRange))))/sigma-iBand*3,[colorArray{count},'--'],'LineWidth',1 )
%     hold on
%     plot(1:length(tempPlot),prcBound(iBand).*ones(length(tempPlot))-iBand*3,[colorArray{count},'-.'] ,'LineWidth',1)
%     hold on
%     plot(1:length(tempPlot),zeros(length(tempPlot))-iBand*3,'k' ,'LineWidth',1)
%     hold on
%     count = count+1 ;
% end
% xlabel('time(ms)')
% set(gca,'ytick',[])
% set(gcf,'Position',[680 127 530 834])

figure;
colorArray = [{'g'},{'r'},{'c'},{'m'},{'k'},{'k'}] ;
count = 1 ;
timeRange = fix(70*fsTemporal)+200:10:fix(74*fsTemporal)+200 ;
plot(1:length(timeRange),zscore(squeeze(sigSmooth(5,5,timeRange)))-3*3,'b','LineWidth',1)
hold on
% plot(1:length(timeRange),zeros(length(tempPlot))-3*3,'k' ,'LineWidth',1)
% hold on

for iBand = 5:2:13
    [temp,mu,sigma] = zscore(real(wvcfs{iBand}(5,5,:))) ; 
    prcPara = 95 ;
    prcBound(iBand) = prctile(abs(temp),prcPara) ;
    tempPlot = temp(timeRange) ;
    plot(1:length(tempPlot),(squeeze(real(tempPlot))) -iBand*3,[colorArray{count},'-'],'LineWidth',1)
    hold on
    plot(1:length(tempPlot),(abs(squeeze(wvcfs{iBand}(5,5,timeRange)))/sigma)-iBand*3,[colorArray{count},'--'],'LineWidth',1 )
    hold on
    plot(1:length(tempPlot),prcBound(iBand).*ones(length(tempPlot))-iBand*3,[colorArray{count},'-.'] ,'LineWidth',1)
    hold on
    plot(1:length(tempPlot),zeros(length(tempPlot))-iBand*3,'k' ,'LineWidth',1)
    hold on
    count = count+1 ;
end
xlabel('time(ms)')
set(gca,'ytick',[])
set(gcf,'Position',[680 127 530 834])

%% fourier domain
figure;
numBand = length(wvcfs) ;
for iBand = 1:numBand
    temp = real(wvcfs{iBand}(5,5,:)) ;
    plot(linspace(0,fsTemporal,length(temp)),abs(fft(squeeze(temp))))
    hold on
end

%% subband frequency and time domain
sigReshape = reshape(sigOri,100,[]) ;
sigReshape2 =  reshape(sigSmooth,400,[]) ;
nfft = 2000 ;
fourCoef = nan(size(sigReshape,1),nfft) ;
fourCoef2 = nan(size(sigReshape2,1),nfft) ;

for iChan = 1:setdiff(1:100,badChannels)
    fourCoef(iChan,:) = abs(fft(sigReshape(iChan,:),nfft)) ;

end
meanFourCoef = nanmean(fourCoef) ;
figure;
loglog(linspace(0,fsTemporal,nfft),meanFourCoef)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([1 300])

for iChan = 190%1:400
    fourCoef2(iChan,:) = abs(fft(sigReshape2(iChan,:),nfft)) ;

end
meanFourCoef2 = nanmean(fourCoef2) ;
figure;
loglog(linspace(0,fsTemporal,nfft),meanFourCoef2,'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([1 300])

%% subband frequency (wavelet)
bandFour = [] ;
for iBand = 1:2:17
bandFour(iBand,:) = abs(fft(real(wvcfs{iBand}(10,10,:)),nfft))/2/pi ;

hold on;
% loglog(linspace(0,fsTemporal,nfft),bandFour(iBand,:),'linewidth',2)
loglog(linspace(0,fsTemporal,nfft),bandFour(iBand,:),'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
end
% xlim([1 200])

%% subband frequency (bandpass)

% bandFour = [] ;
% for iBand = 1:2:17
% bandFour(iBand,:) = abs(fft((bandpassSig(iBand,10,10,:)),nfft))/2/pi ;
% % loglog(linspace(0,fsTemporal,nfft),bandFour(iBand,:),'linewidth',2)
% loglog(linspace(0,fsTemporal,nfft),bandFour(iBand,:),'linewidth',2)
% hold on;
% 
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% end

%% using chronus toolbox
S = [] ;
param.pad = 3 ;
param.Fs = fsTemporal ;
% for iChan = 1:setdiff(1:100,badChannels)
%     [S(iChan,:),f] = mtspectrumc(sigReshape(iChan,:),param) ;
% 
% end
% meanFourCoef = nanmean(S) ;

[S,f] = mtspectrumc(sigReshape2(190,:),param) ;
figure;
loglog(f,S,'g','lineWidth',1)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([1 300])

for iBand = 1:2:17
    S = [] ;
    param.pad = 3 ;
    param.Fs = fsTemporal ;
    [S,f] = mtspectrumc(real(wvcfs{iBand}(:)),param) ;

    hold on;
    % meanFourCoef = nanmean(S) ;
    loglog(f,S,'.','lineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
end
%%

%%
clearvars LFPs sigOri sigSmooth wvcfs temp