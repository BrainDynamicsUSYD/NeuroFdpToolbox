dataFileName = 'ma027_032' ;
timeStart = 0 ;     % for finding centres
timeStart2 = 0;
subBand = [30 , 80] ;
plotAmp = 1 ;
% timeStep = 5 ;
% totalSur = 200 ;
minBurstTime = 30 ;
% minSize = 4;

load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

sigOriTemp = reshape(sigOri,10*10,[]) ;
surSigTemp = nan(size(sigOriTemp)) ;

randNumHalf =  2*pi*rand(size(sigOri,3)/2-1,1) ;
randNum = [0;randNumHalf;0;-flip(randNumHalf) ]' ;
for iChannel = setdiff(1:100,badChannels)
    freqSig = fft(sigOriTemp(iChannel,:)) ;
    absFreq = abs(freqSig) ;
    phaseFreq = angle(freqSig) ;
    reconSig = ifft(absFreq.*exp(1i*phaseFreq)) ;
    
    surSigTemp(iChannel,:) = ifft(absFreq.*exp(1i*(phaseFreq+randNum))) ;
end
surSig = reshape(real(surSigTemp),10,10,[]) ;


% find subband signals
bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels,0) ;
% bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;

if plotAmp
    sigIn = abs(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:100*fsTemporal) ));
else
    sigIn = angle(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:100*fsTemporal) ));
end

%%
sigIn = abs(squeeze(hilbertSig(1,5,5,0*fsTemporal+1:100*fsTemporal) ));

%% 95 prctile
boundPrctile = 95 ;
prcBound = prctile(sigIn(:),boundPrctile) ;
sigBinary = zeros(size(sigIn)) ;
sigBinary(sigIn>prcBound) = 1;
sigPlot = sigBinary.*sigIn ;

%% Yifan
sigInReshape = reshape(sigIn,100,[]) ;
gaus_width = 12.5; %  ms;
[Kernel] = spike_train_kernel_YG(gaus_width,1/fsTemporal*1e3,'gaussian_unit');
for iChannel = 1:100
        sigIntemp(iChannel,:) = conv(sigInReshape(iChannel,:),Kernel,'same');
end
GammaBurstEvent = find_Burst_1D(sigIntemp,fsTemporal,0,badChannels,100) ; 

% GammaBurstEvent2 = find_Burst_1D(sigInReshape,fsTemporal,0,badChannels,100) ; 


%% duration distribution for Yifan
sigInProp = [] ;
% sigInProp = GammaBurstEvent2.burst_du_steps{:} ;
for iChannel = 1:99
    sigInProp = [sigInProp,GammaBurstEvent.burst_du_steps{iChannel}] ;  % duration
    % sigInProp = [sigInProp, GammaBurstEvent.flat_du_steps{iChannel}] ;
end

%% lognoraml fit
close all
histogram(sigInProp,100)
title('Duration distribution of single electrode')
xlabel('Time (ms)')
figure;

[PARMHAT,b] =  lognfit(sigInProp) ;

[n,x] = histcounts(sigInProp,200,'normalization','pdf') ;
semilogx((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = lognpdf(xNew,PARMHAT(1),PARMHAT(2)) ;
semilogx(xNew,y,'r')
title('Duration distribution of single electrode')
xlabel('Time (ms)')
ylabel('probability')

legend('original data','log-normal fit')

% alpha stable fit
figure
pd = fitdist(sigInProp','stable') ;

[n,x] = histcounts(sigInProp,200,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
semilogy(xNew,y,'r')
title('Duration distribution of single electrode')
xlabel('Time (ms)')
ylabel('probability')

legend('original data','\alpha stable distribution fit')

% power law fit
figure;
[n, xout] = histcounts(sigInProp,200,'normalization','pdf');
% n=n/sum(n) ;
loglog((xout(2:end)+xout(1:end-1))/2,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Duration distribution of single electrode')
xlabel('Time (ms)')
ylabel('Probability')
nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(2:1*length(nFit))),log(nFit(2:1*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
% xlim([20 3e3])   
% ylim([1e-5 1])
legend('original data','power law fitting')

% exponential fit
figure;
pd = fitdist(sigInProp','exp') ;

[n,x] = histcounts(sigInProp,200,'normalization','pdf') ;
plot((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
plot(xNew,y,'r')
title('Duration distribution of single electrode')
xlabel('Time (ms)')
ylabel('probability')

legend('original data','exponential distribution fit')

%% duration distribution for 95 prctile
duration = [] ;
interVal = [] ;
for iChannel = setdiff(1:100,badChannels)
    CC = bwconncomp(squeeze(sigBinary(mod(iChannel-1,10)+1,ceil(iChannel/10) ,:))) ;                % 3D
    B1 = regionprops(CC,'BoundingBox');
    boundary1 = cat(1, B1.BoundingBox);
% distribution of duration
    duration = [duration;boundary1(2:end-1,4)] ;

% distribution of interval
    interVal = [interVal;(boundary1(3:end-1,2)-boundary1(2:end-2,2)-boundary1(2:end-2,4))] ;
end

%%
close all
hist(duration,40)
figure;
hist(interVal,200)


%% burst detection (15ms/30ms threshold)

%% wavelet transform
%%
close all
for timeStep = 25
    addpath('ToolOthers/uimage')
    timeEpoch = 2 ;
    timeStart = timeEpoch*(timeStep-1)  ;
    timeEnd = timeEpoch*timeStep   ;
    %timeStart = 0 ;
    %timeEnd = 20 ;
    
    freqRange = 20:80 ;
    waveletParam = 'cmor1.5-1' ;
    fc = centfrq(waveletParam) ;
    scalerange = fc./(freqRange/fsTemporal) ;
    scales = scalerange(end):0.5:scalerange(1) ;
    pseudoFreq = scal2frq(scales, waveletParam, 1/fsTemporal) ;
    
    sigRange = floor(timeStart*fsTemporal)+1 : floor(timeEnd*fsTemporal) ;
    tempData =  squeeze (sigIn(sigRange)) ; 
%     tempData = squeeze(dataFull( 8,floor(timeStart*fsTemporal)+1 : ...
%         floor(timeEnd*fsTemporal) ))  ;
    
    wt = cwt( tempData ,scales, waveletParam  ) ;
    tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
    
%     figure
%     wtNorm = zscore(abs(wt),[],2) ;
%     
%     imagesc(tempTimeAxis,pseudoFreq,(wtNorm(:,:)))
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
    subplot(4,1,1)
    uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(abs(wt(end:-1:1,:)) ))
    % uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),(zscore(abs(wt(end:-1:1,:)) )))
    % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    set(gca,'YDir','normal')
    title('Wavelet spectrum')
end

%
subplot(4,1,2)
plot(tempTimeAxis,GammaBurstEvent.is_burst(sigRange))
xlabel('Time (s)')
title('Smoothed 2SD burst detection')
%
subplot(4,1,3)
plot(tempTimeAxis,GammaBurstEvent2.is_burst(sigRange))
xlabel('Time (s)')
title('2SD Burst detection')

%
subplot(4,1,4)
plot(tempTimeAxis,sigIn(sigRange))
hold on
plot(tempTimeAxis,sigIntemp(sigRange),'r--')
xlabel('Time (s)')
legend('Hilbert envelope', 'Smoothed envelope')

%%
CentroidsFull = [];
for xChannel = 2:9
    for yChannel = 2:9
for timeStep = 1
    % addpath('ToolOthers/uimage')
    timeEpoch = 100 ;
    timeStart = timeEpoch*(timeStep-1)  ;
    timeEnd = timeEpoch*timeStep   ;
    %timeStart = 0 ;
    %timeEnd = 20 ;
    
    freqRange = 30:80 ;
    waveletParam = 'cmor1.5-1' ;
    fc = centfrq(waveletParam) ;
    scalerange = fc./(freqRange/fsTemporal) ;
    scales = scalerange(end):0.5:scalerange(1) ;
    pseudoFreq = scal2frq(scales, waveletParam, 1/fsTemporal) ;
    
    sigRange = floor(timeStart*fsTemporal)+1 : floor(timeEnd*fsTemporal) ;
    tempData =  squeeze (sigIn(xChannel,yChannel,sigRange)) ; 
%     tempData = squeeze(dataFull( 8,floor(timeStart*fsTemporal)+1 : ...
%         floor(timeEnd*fsTemporal) ))  ;
    
    wt = cwt( tempData ,scales, waveletParam  ) ;
    tempTimeAxis = linspace(timeStart,timeEnd,size(wt,2)) ;
end

%%CentroidsFull
sigWTabs = abs(wt) ;
prc = prctile(sigWTabs,95) ;
binaryWT = zeros(size(sigWTabs)) ;
binaryWT(sigWTabs>prc) = 1 ;
CC = bwconncomp(binaryWT) ;                % 3D
B1 = regionprops(CC,'BoundingBox');

S = regionprops(CC,'Centroid');
Centroids = cat(1, S.Centroid);
CentroidsFull = [CentroidsFull;Centroids] ;
    end
end


%%
cFreq = pseudoFreq(round(CentroidsFull(:,2))) ;
histogram(cFreq,25)
