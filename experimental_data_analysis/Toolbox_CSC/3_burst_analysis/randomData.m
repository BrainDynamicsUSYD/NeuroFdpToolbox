clc 
clear

%% random data test
fsTemporal = 1000;
x = (rand(10000,1)) ;

freqRange = [1,80] ;
scaleStep = 0.5 ;
wName = 'cmor1.5-1' ;
[cwtCoef,pseudoFreq] = find_cwtCoef2(x', freqRange, scaleStep,...
   fsTemporal,wName) ;


% zscoreCoef =abs(cwtCoef) ;
figure
% uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
uimagesc(1:10000,pseudoFreq(end:-1:1),  abs(cwtCoef(end:-1:1,:)) )
% 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(gca,'YDir','normal')

%%

x = rand(100000,1) ;

freqRange = [30,80] ;
scaleStep = 0.5 ;
wName = 'cmor1.5-1' ;
[cwtCoef,pseudoFreq] = find_cwtCoef2(x', freqRange, scaleStep,...
   fsTemporal,wName) ;

absCoef =abs(cwtCoef) ;
figure
% uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
uimagesc(1:length(x),pseudoFreq(end:-1:1),  absCoef(end:-1:1,:) )
% 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(gca,'YDir','normal')


%%

Y = prctile(absCoef(:),95) ;
binaryImage = zeros(size(absCoef)) ;
binaryImage(absCoef>=Y) = 1 ;
%
CC = bwconncomp(binaryImage(:,:)) ;
validRegionIdx = 0;
count = 1;
for iRegion = 1:size(CC.PixelIdxList,2)
    if(size(CC.PixelIdxList{iRegion},1)>200)
        validRegionIdx(count) = iRegion ;
        count = count + 1 ;
    end
end

S = regionprops(CC,'Centroid');
centroids = cat(1, S.Centroid);


% cTime = centroids(validRegionIdx(:),1)/fsTemporal ;
% cFreq = pseudoFreq(round(centroids(validRegionIdx(:),2))) ;
% 
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

hist(boundary(validRegionIdx(:),3),20)


%%
figure
% uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
uimagesc(1:10000,pseudoFreq(end:-1:1),  binaryImage(end:-1:1,:) )
% 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(gca,'YDir','normal')


%% bandpass signals

subBand = [30, 80] ;
bandpassSig = find_bandpassSig(x',subBand, fsTemporal,2, [],0,8,0) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4,[],0) ;

absHilbert = squeeze(abs(hilbertSig)) ;
GammaBurstEvent = find_Burst_1D(absHilbert',fsTemporal) ;

%%
hist(GammaBurstEvent.burst_du_steps{:},20)

%% spatial pattern
xMatrix = rand(10,10,100000) ;
fsTemporal = 1000 ;

% find subband signals
bandpassSig = find_bandpassSig(xMatrix,subBand, fsTemporal,3) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;

plotAmp = 1 ;
if plotAmp
    sigIn = abs(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:90*fsTemporal) ));
else
    sigIn = angle(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:90*fsTemporal) ));
end

clearvars bandpassSig hilbertSig


%% Interpolation and smoothing

close all
plotLength = fix(10*fsTemporal) ;
sigSmooth = zeros(20,20,plotLength) ;
timeStart = 5000 ;

addpath(genpath([pwd,'/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))

resizeScale = 2 ;
for iTime = 1:plotLength
    timeSlot = iTime+timeStart ;
     %Try smoothing
    filtWidth = 3;
    filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);   
end
clearvars smoothTemp
%% find 95 percentile
if ~plotAmp
    sigBinary = zeros(size(sigSmooth)) ;
    sigBinary(sigSmooth>pi*5/6) = 1 ;
    sigPlot = sigBinary.*sigSmooth ;
else
    prc95 = prctile(sigSmooth(:),95) ;
    sigBinary = zeros(size(sigSmooth)) ;
    sigBinary(sigSmooth>prc95) = 1;
    sigPlot = sigBinary.*sigSmooth ;
end

%% find central
largestCentroids = nan(size(sigBinary,3),2) ;
boundaryRaw = nan(size(sigBinary,3),4) ;
% count = 1;
for iTime = 1:size(sigBinary,3)
    CC = bwconncomp(sigBinary(:,:,iTime)) ;

    if (size(CC.PixelIdxList,2) == 0)
        continue
    else
        sizeRegion = zeros(size(CC.PixelIdxList,2),1) ;
        for iRegion = 1:size(CC.PixelIdxList,2)
            sizeRegion(iRegion) = size(CC.PixelIdxList{iRegion},1) ;
        end
        [regionSize,largestRegionIdx] = max(sizeRegion) ;
        if regionSize<4
            continue
        else
            S = regionprops(CC,'Centroid');
            Centroids = cat(1, S.Centroid);
            
            largestCentroids(iTime,:) = Centroids(largestRegionIdx,:) ;
            % count = count+1 ;
            B = regionprops(CC,'BoundingBox');
            boundaryTemp = cat(1, B.BoundingBox);
            boundaryRaw(iTime,:) = boundaryTemp(largestRegionIdx,:) ;
        end
    end
end

% %% calculate the duration of each cluster
% plot(largestCentroids(:,1))
% hold on
% % plot(largestCentroids(:,2))-20
% % first find drift location
% diffCent = diff(largestCentroids,1) ;
% driftLoc = find((abs(diffCent(:,1))>10) | (abs(diffCent(:,2))>10) )  ;
% 
% % find continuous bursts
% binaryCentTemp = zeros(size(largestCentroids,1),1) ;
% binaryCentTemp(~isnan(largestCentroids(:,1))) = 1 ;
% binaryCent = binaryCentTemp ;
% 
% for iDrift = 1:length(driftLoc)
%     binaryCent = [binaryCent(1:(driftLoc(iDrift)+iDrift-1));0;...
%         binaryCent((driftLoc(iDrift)+iDrift):end)] ;
% end
% plot(binaryCent)

%%
CC = bwconncomp(binaryCent) ;
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);
hist(boundary(:,4),20)

%% Centroids with long duration
% first find drift location
diffCent = diff(largestCentroids,1) ;
driftLoc = find((abs(diffCent(:,1))>4) | (abs(diffCent(:,2))>4) )  ;
% find continuous bursts
binaryCentTemp2 = zeros(size(largestCentroids,1),1) ;
binaryCentTemp2(~isnan(largestCentroids(:,1))) = 1 ;
binaryCentTemp2(driftLoc) = 0 ;

CC = bwconncomp(binaryCentTemp2) ;
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

minDuration = 30/1000*fsTemporal ;
selectCentroids = nan(size(largestCentroids)) ;
selectBoundary = nan(size(boundaryRaw)) ;
validIdx = find(boundary(:,4)>minDuration) ;
for idx = 1:length(validIdx)
    selectCentroids(CC.PixelIdxList{validIdx(idx)},:) = 1 ;
    selectBoundary(CC.PixelIdxList{validIdx(idx)},:) = 1 ;
end
%
largestCentroids = largestCentroids.*selectCentroids ;
boundaryValid = boundaryRaw.*selectBoundary;
%%
areaBurst = boundaryValid(:,3).*boundaryValid(:,4)/4 ;
areaBurst(isnan(areaBurst)) = [] ;
meanArea = num2str(mean(areaBurst)) ;

hist(areaBurst,20)
title(['Histogram of burst area, mean = ', meanArea])
xlabel('Area (electrode^2)')
%% Calculate the distance of Gamma burst
% electrode distance 0.4mm, 0.2mm after interpolation 
distCent = 0.5*sqrt(sum((largestCentroids(2:end,:) - largestCentroids(1:end-1,:)).^2,2)) ;
% distCent = (largestCentroids(2:end,2) - largestCentroids(1:end-1,2)) ;
distCentMM = 0.4*distCent ;
% hist(distanceCent,40)
% xlabel('Distance (mm)')

%% Histogram of Gamma drift distsance
close all
[n, xout] = hist(distCent,200);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% xlabel('Distance (unit electrode)')
xlabel('Distance (electrode)')
%xlim([0.1,15])
% xlim([0,14*0.4])
hold on

% % fit a levy distribution
% pd = fitdist(distCent,'stable') ;
% x_values = -10:1:10;
% y = pdf(pd,x_values);
% plot(x_values,y,'r','LineWidth',2)


%% make a movie
makeMovie = 1;
timeStart2 = 0 ;
if makeMovie
    vidTitle = [pwd,'/Results/GammaBurst/','Full','randomNoise',num2str(timeStart2)] ;
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = 20 ;
    open(vidObj);
     fig=figure ;
     set(gcf,'Position',[260 23 1159 926])
%     
    
    for iTime = 1:1000
        timeSlot = iTime+timeStart2 ;
        subplot(2,2,1)
        imagesc(sigIn(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        title(['Raw signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.8*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,2)
        imagesc(sigSmooth(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Smoothed signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.8*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,3)
        imagesc(sigPlot(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Detected burst at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.8*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,4)
        imagesc(sigSmooth(:,:,iTime))
        currentPoint = largestCentroids(iTime,:) ;
        % nextPoint = largestCentroids(iTime+1,:) ;
        plot(currentPoint(:,1),currentPoint(:,2),'.','MarkerSize',20)
        xlim([0,20])
        ylim([0,20])
        title(['Detected burst Centre at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        colorbar
        
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        % pause(0.01)
        
    end
    close(vidObj);
end