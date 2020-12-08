%% time hierachy
fsTemporal = 1017.25262451172;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

% find subband signals
subBand = [30,80] ;
bandpassSig = find_bandpassSig(sigOri,subBand, fsTemporal,3) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;

sigIn = abs(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:150*fsTemporal) ));

%%
% close all

sigOriPick = [] ;
n1 = 5;
for iRep = 1:n1
    sigOriPick = [sigOriPick;squeeze(sigIn(5,5,iRep:n1:end))] ;
end
diffSig = diff(sigOriPick) ;
figure
[x,n] = hist(diffSig,400);
plot(n,x/sum(x)/10^(1.5),'go');

sigOriPick = [] ;
n2 = 20;
for iRep = 1:n2
    sigOriPick = [sigOriPick;squeeze(sigIn(5,5,iRep:n2:end))] ;
end
diffSig = diff(sigOriPick) ;
hold on
[x,n] = hist(diffSig,400);
plot(n,x/sum(x)/10^1,'bs');

sigOriPick = [] ;
n3 = 80 ;
for iRep = 1:n3
    sigOriPick = [sigOriPick;squeeze(sigIn(5,5,iRep:n3:end))] ;
end
diffSig = diff(sigOriPick) ;
hold on
[x,n] = hist(diffSig,400);
plot(n,x/sum(x)/10^(0.5),'cd');

sigOriPick = [] ;
n4= 400 ;
for iRep = 1:n4
    sigOriPick = [sigOriPick;squeeze(sigIn(5,5,iRep:n4:end))] ;
end
diffSig = diff(sigOriPick) ;
hold on
[x,n] = hist(diffSig,400);
plot(n,x/sum(x),'rp');

legend([num2str(n1),'ms interval'],[num2str(n2),'ms interval'],...
    [num2str(n3),'ms interval'],[num2str(n4),'ms interval'])
title('Raw signal variance in one electrode time hierachy')
xlabel('amplitude variance')
ylabel('probability of amplitude variance')
set(gca, 'YScale', 'log')
%xlim([-5e-5 5e-5])

%%

sigOriPick = squeeze(sigOri(5,5,:)) ;
diffSig = diff(sigOriPick) ;
figure
[x,n3] = hist(diffSig,400);
plot(n3,x,'go');
hold on

% fit a levy distribution
pd = fitdist(diffSig(diffSig>-1.5e-5 & diffSig<1.5e-5),'stable') ;
x_values = -5e-5:0.01e-5:5e-5;
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',2)

title('instantaneous amplitude variance time hierachy')
xlabel('amplitude variance')
ylabel('probability of amplitude variance')
set(gca, 'YScale', 'log')

%% 
sigOriPick = [] ;
n3= 125 ;
for iRep = 1:n3
    sigOriPick = [sigOriPick;squeeze(sigOri(5,5,iRep:n3:end))] ;
end
diffSig = diff(sigOriPick) ;
hold on
[x,n3] = hist(diffSig,400);
plot(n3,x/sum(x),'bp');

% fit a levy distribution
pd = fitdist(diffSig(diffSig>-2e-5 & diffSig<2e-5),'stable') ;
x_values = -3e-4:0.01e-4:3e-4;
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',2)

title('instantaneous amplitude variance time hierachy')
xlabel('amplitude variance')
ylabel('probability of amplitude variance')
set(gca, 'YScale', 'log')


%% autocorrelation
close all
autocorr(patternScale)
title('Autocorrelation\_Scale\_my144\_30ms')
xlabel('Lag (burst number)')

figure
autocorr(Duration)
title('Autocorrelation\_Duration\_my144\_30ms')
xlabel('Lag (burst number)')
%%
figure
temp = instantScale(1:5:end,3) ;
temp = temp(temp~=0) ;
autocorr(temp)
title('Autocorrelation\_SitesInBurst\_my144\_30ms')
xlabel('Lag (5 ms)')

%% cross correlation
scaleInBurst = [] ;
nPoints = 100 ;
for iBurst = 1:size(instantScale,2)
    temp = instantScale(:,iBurst) ;
    temp = temp(temp~=0) ;
    if length(temp) <nPoints
        temp = [temp;zeros(nPoints-length(temp)+1,1)] ;
    end
    scaleInBurst = [scaleInBurst,temp(1:nPoints)] ;
end
r = corrcoef(scaleInBurst) ;
rPlot = r ;
rPlot(rPlot<0.9) = 0 ;
imagesc(rPlot)

%% sticky and flights

for iBurst = 1:40   % 24
    displacementC = [] ;
    % center = squeeze(Centroids(:,iBurst,:)) ;
    % center = center(center(:,1)~=0,:) ;
    center = WCentroids{iBurst} ;
    for iTime = 1:size(center,1)
        displacementC(iTime) = sum((center(1,:)-center(iTime,:)).^2) ;
    end
    plot(displacementC)
    pause
    close all
end

%%
numFlights = 0 ;
numSticky = 0 ;
flightsDuration = [] ;
stickyDuration = [] ;

for iBurst = 1:size(Centroids,2) 
    displacementC = [] ;
    % center = squeeze(WCentroids(:,iBurst,:)) ;
    % center = center(center(:,1)~=0,:) ;
    center = WCentroids{iBurst} ;
    for iTime = 1:size(center,1)
        displacementC(iTime) = sum((center(1,:)-center(iTime,:)).^2) ;
    end
    [PKS,LOC] = findpeaks(displacementC) ;
    diffPKS = diff(PKS) ;
    idxFlights = (abs(diffPKS)<=1.5) ;
    diffIdx = diff(idxFlights) ;
    flightsStart = find(diffIdx==1) ;
    flightsEnd = find(diffIdx==-1) ;
    if(isempty(flightsStart))
        %flightsDuration{iBurst} = flightsEnd - 1 ;
        continue
    end
    if(isempty(flightsEnd))
        %flightsDuration{iBurst} = length(displacementC) - flightsStart ;
        continue
    end
    
    if (flightsStart(1))>(flightsEnd(1))
        flightsEnd(1) = [] ;
    end
    
    if length(flightsStart)<length(flightsEnd)
        error('error!') ;
    elseif length(flightsStart)>length(flightsEnd)
        flightsDuration{iBurst} = flightsEnd - flightsStart(1:end-1) ;
    else
        flightsDuration{iBurst} = flightsEnd - flightsStart ;
    end
    
    idxSticky = (abs(diffPKS)>1.5) ;
    diffIdx = diff(idxSticky) ;
    stickyStart = find(diffIdx==1) ;
    stickyEnd = find(diffIdx==-1) ;
    if (stickyStart(1))>(stickyEnd(1))
        stickyEnd(1) = [] ;
    end
    if length(stickyStart)<length(stickyEnd)
        error('error!') ;
    elseif length(stickyStart)>length(stickyEnd)
        stickyDuration{iBurst} = stickyEnd - stickyStart(1:end-1) ;
    else
        stickyDuration{iBurst} = stickyEnd - stickyStart ;
    end
end
figure
temp = cell2mat(flightsDuration) ;
[x,n]=hist(temp,40) ;
loglog(n,x,'o')
xlabel('t(ms)')
ylabel('count')
title('flights time distribution of my144')

figure
temp2 = cell2mat(stickyDuration) ;
[x,n]=hist(temp2,40) ;
loglog(n,x,'o')
xlabel('t(ms)')
ylabel('count')
title('sticking time distribution of my144')

%% analysis within burst
% distribution of scale

% scale and distance correlation
% total power/ peak power and distance correlation
center = WCentroids ;
for iBurst = 1:size(Centroids,2)
    scaleBefFlig = instantScale(:,iBurst) ;
    scaleBefFlig = scaleBefFlig(scaleBefFlig~=0) ;

    centreBurst = center(center(:,iBurst,1)~=0,iBurst,:) ;
    centreBurst = squeeze(centreBurst) ;
    % distCentBurst = sqrt( sum( ( centreBurst(1:end-1,:)-centreBurst(2:end,:) ).^2,2) ) ;
    distCentBurst = sqrt( sum( (diff(centreBurst)).^2,2 ) ) ;
    
    [r,p] = corrcoef(scaleBefFlig(1:end-1),distCentBurst) ;
    [r2,p2] = corrcoef(sumAmp(1:end-1),distCentBurst) ;
    [r3,p3] = corrcoef(peakAmp(1:end-1),distCentBurst) ;
    
    dist_scale_cor(iBurst) = r(2) ;
    dist_SumPower_cor(iBurst) = r2(2) ;
    dist_PeakAmp_cor(iBurst) = r3(2) ;
end
plot(dist_scale_cor)
figure
plot(dist_SumPower_cor)
figure
plot(dist_PeakAmp_cor)

% scale and duration correlation

% power and duration correlation

%%
numFlights = 0 ;
numSticky = 0 ;
flightsDuration = [] ;
stickyDuration = [] ;
scaleBefFlig = [] ;
distFlight = [] ;

for iBurst = 1:size(Centroids,2) 
    displacementC = [] ;
    center = squeeze(WCentroids(:,iBurst,:)) ;
    center = center(center(:,1)~=0,:) ;
    for iTime = 1:size(center,1)
        displacementC(iTime) = sum((center(1,:)-center(iTime,:)).^2) ;
    end
    [PKS,LOC] = findpeaks(displacementC) ;
    diffPKS = diff(PKS) ;
    idxFlights = (abs(diffPKS)<=4) ;
    diffIdx = diff(idxFlights) ;
    flightsStart = find(diffIdx==1) ;
    flightsEnd = find(diffIdx==-1) ;
    if(isempty(flightsStart))
        flightsDuration{iBurst} = flightsEnd - 1 ;
        continue
    end
    if(isempty(flightsEnd))
        flightsDuration{iBurst} = length(displacementC) - flightsStart ;
        continue
    end
    
    if (flightsStart(1))>(flightsEnd(1))
        flightsEnd(1) = [] ;
    end
    if (flightsStart(1)==1)
        flightsStart(1) = [] ;
        flightsEnd(1) = [] ;
    end
    instantScaleBurst = instantScale(instantScale(:,iBurst)~=0,iBurst) ;

    if length(flightsStart)<length(flightsEnd)
        error('error!') ;
    elseif length(flightsStart)>length(flightsEnd)
        startPoint = flightsStart(1:end-1) ;
        flightsDuration{iBurst} = flightsEnd - startPoint ;
        scaleBefFlig{iBurst} = instantScaleBurst(startPoint-1) ;
        distFlight{iBurst} = sum( (center(startPoint,:)...
            -center(startPoint-1,:)) .^2) ;
    else
        flightsDuration{iBurst} = flightsEnd - flightsStart ;
        scaleBefFlig{iBurst} = instantScaleBurst(flightsStart-1) ;
        distFlight{iBurst} = sum( (center(flightsStart,:)...
            -center(flightsStart-1,:)) .^2) ;
        
    end
    
    
    idxSticky = (abs(diffPKS)>4) ;
    diffIdx = diff(idxSticky) ;
    stickyStart = find(diffIdx==1) ;
    stickyEnd = find(diffIdx==-1) ;
    if (stickyStart(1))>(stickyEnd(1))
        stickyEnd(1) = [] ;
    end
    if length(stickyStart)<length(stickyEnd)
        error('error!') ;
    elseif length(stickyStart)>length(stickyEnd)
        stickyDuration{iBurst} = stickyEnd - stickyStart(1:end-1) ;
    else
        stickyDuration{iBurst} = stickyEnd - stickyStart ;
    end
    
end

[r,p] = corrcoef(scaleBefFlig{1},flightsDuration{1}) ;

[r2,p2] = corrcoef(scaleBefFlig{1},distFlight{1}) ;


%% plot bursts
center = WCentroids ;
for iBurst = 8% 1:size(Centroids,2)
    scaleBefFlig = instantScale(:,iBurst) ;
    scaleBefFlig = scaleBefFlig(scaleBefFlig~=0) ;

    centreBurst = center(center(:,iBurst,1)~=0,iBurst,:) ;
    centreBurst = squeeze(centreBurst) ;
    
     plot(centreBurst(:,1),centreBurst(:,2),'o-')
     xlim([0 20])
     ylim([0 20])
     hold on
end