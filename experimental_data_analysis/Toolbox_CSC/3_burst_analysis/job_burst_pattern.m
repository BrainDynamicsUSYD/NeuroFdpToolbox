% OLD method
function job_burst_pattern(arrayID)

cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%
% arrayID = 8 ;

switch arrayID
    case 1
        dataFileName = 'my144_101' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;        % for movie
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;      
        
    case 2
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 1;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 3
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 2;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 1 ;
        
    case 4
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;        % for movie
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;
        
    case 5
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 1;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 6
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 1 ;
        
    case 7
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;
        
    case 8
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 1;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 9
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 1 ;
        
    case 10
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;
        
    case 11
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 1;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 12
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 1 ;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
clearvars LFPs
sigOriTemp = reshape(sigOri,10*10,[]) ;
surSigTemp = nan(size(sigOriTemp)) ;

%%
totalSur = 1 ;
for nSur = 1:totalSur
%     randNumHalf =  2*pi*rand(size(sigOri,3)/2-1,1) ;
%     randNum = [0;randNumHalf;0;-flip(randNumHalf) ]' ;
%     for iChannel = setdiff(1:100,badChannels)
%         freqSig = fft(sigOriTemp(iChannel,:)) ;
%         absFreq = abs(freqSig) ;
%         phaseFreq = angle(freqSig) ;
%         reconSig = ifft(absFreq.*exp(1i*phaseFreq)) ;
%         
%         surSigTemp(iChannel,:) = ifft(absFreq.*exp(1i*(phaseFreq+randNum))) ;
%     end
%  surSig = reshape(real(surSigTemp),10,10,[]) ;

surSig = sigOri ;

% find subband signals
bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;

if plotAmp
    sigIn = abs(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:150*fsTemporal) ));
else
    sigIn = angle(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:150*fsTemporal) ));
end

clearvars bandpassSig hilbertSig


% %% Interpolation and smoothing
% 
% close all
% plotLength = fix(100*fsTemporal) ;
% sigSmooth = zeros(20,20,plotLength) ;
% % timeStart = 5000 ;
% 
% addpath(genpath([pwd,'/ToolNeuroPatt']))
% addpath(genpath([pwd,'/ToolOthers/nanconv']))
% 
% resizeScale = 2 ;
% for iTime = 1:plotLength
%     timeSlot = iTime+timeStart ;
%      %Try smoothing
%     filtWidth = 3;
%     filtSigma = 1;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
%     sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);   
% end
% clearvars smoothTemp


%% find 95 percentile
sigSmoothTemp = sigIn ;
if ~plotAmp
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>pi*5/6) = 1 ;
    sigPlot = sigBinary.*sigSmoothTemp ;
else
    boundPrctile = 90 ;
    prcBound = prctile(sigSmoothTemp(:),boundPrctile) ;
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>prcBound) = 1;
    sigPlot = sigBinary.*sigSmoothTemp ;
end

%% find central
largestCentroids = nan(size(sigBinary,3),2) ;
boundaryRaw = nan(size(sigBinary,3),4) ;
% count = 1;

% loop through all the time frames to find patterns
for iTime = 1:size(sigBinary,3)
    % find connected points (2D) in each time frame 
    CC = bwconncomp(sigBinary(:,:,iTime)) ;

    if (size(CC.PixelIdxList,2) == 0)
        continue
    else
        sizeRegion = zeros(size(CC.PixelIdxList,2),1) ;
        for iRegion = 1:size(CC.PixelIdxList,2)
            sizeRegion(iRegion) = size(CC.PixelIdxList{iRegion},1) ;
        end
        [regionSize,largestRegionIdx] = max(sizeRegion) ;
        if regionSize<minSize
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
% CC = bwconncomp(binaryCent) ;
% B = regionprops(CC,'BoundingBox');
% boundary = cat(1, B.BoundingBox);
% hist(boundary(:,4),20)

%% Centroids with long duration (guarantee a full 30Hz Gamma cycles)
% first find drift location
diffCent = diff(largestCentroids,1) ;
driftLoc = find((abs(diffCent(:,1))>4) | (abs(diffCent(:,2))>4) )  ;
% find continuous bursts (get rid of drifts by setting the end point to zero)
binaryCentTemp2 = zeros(size(largestCentroids,1),1) ;
binaryCentTemp2(~isnan(largestCentroids(:,1))) = 1 ;
binaryCentTemp2(driftLoc) = 0 ;

% detect continous time series 
CC = bwconncomp(binaryCentTemp2) ;
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

minDuration = minBurstTime/1000*fsTemporal ;
validCentroids = nan(size(largestCentroids)) ;
validBoundary = nan(size(boundaryRaw)) ;
validIdx = find(boundary(:,4)>minDuration) ;

burstDuration = zeros(length(validIdx),totalSur) ;
burstXwidth = zeros(length(validIdx),totalSur) ; 
burstYwidth = zeros(length(validIdx),totalSur) ;
burstArea = zeros(length(validIdx),totalSur) ;

peak = zeros(length(validIdx),totalSur) ;
timeDist = zeros(length(validIdx) - 1,totalSur) ;
spaceDist = zeros(length(validIdx) - 1,totalSur) ;

for idx = 1:length(validIdx)
    currentIdx = CC.PixelIdxList{validIdx(idx)} ;
    currentBurstB = boundaryRaw(currentIdx,:)   ;
    burstDuration(idx,nSur) = boundary(validIdx(idx),4) ;
    burstXwidth(idx,nSur) = sum(currentBurstB(:,3))/burstDuration(idx,nSur) ;
    burstYwidth(idx,nSur) = sum(currentBurstB(:,4))/burstDuration(idx,nSur) ;
    burstArea(idx,nSur) = burstXwidth(idx,nSur)* burstYwidth(idx,nSur) ;
    % size(currentBurst,1) == burstDuration(idx)
    
    validBoundary(currentIdx,:) = 1 ;
    validCentroids(currentIdx,:) = 1 ; 
    
    currentBurstAmp = (sigSmoothTemp(:,:,currentIdx).*sigBinary(:,:,currentIdx)) ;
    peak(idx,nSur) = max(currentBurstAmp(:)) ;
    if idx<length(validIdx)
        nextIdx = CC.PixelIdxList{validIdx(idx+1)} ;
        timeDist(idx,nSur) = nextIdx(1) - (currentIdx(end)+1) + 1 ;
        spaceDist(idx,nSur) = sqrt(sum((largestCentroids(currentIdx(end,:)) + ...
            largestCentroids(nextIdx(1,:))).^2)) ;
    end   
end
%
% largestCentroids = largestCentroids.*validCentroids ;
% boundaryValid = boundaryRaw.*validBoundary;
%%
% areaBurst = boundaryValid(:,3).*boundaryValid(:,4)/4 ;
% areaBurst(isnan(areaBurst)) = [] ;
% meanArea(nSurr) = (mean(areaBurst)) ;




% %% Calculate the distance of Gamma burst
% % electrode distance 0.4mm, 0.2mm after interpolation 
% distCent = 0.5*sqrt(sum((largestCentroids(2:end,:) - largestCentroids(1:end-1,:)).^2,2)) ;
% % distCent = (largestCentroids(2:end,2) - largestCentroids(1:end-1,2)) ;
% distCentMM = 0.4*distCent ;
% % hist(distanceCent,40)
% % xlabel('Distance (mm)')
% 
% %% Histogram of Gamma drift distsance
% close all
% [n, xout] = hist(distCent,200);
% bar(xout, n, 'barwidth', 1, 'basevalue', 1);
% % set(gca,'XScale','log')
% % set(gca,'YScale','log')
% % xlabel('Distance (unit electrode)')
% xlabel('Distance (electrode)')
% %xlim([0.1,15])
% % xlim([0,14*0.4])
% hold on
% 
% % % fit a levy distribution
% % pd = fitdist(distCent,'stable') ;
% % x_values = -10:1:10;
% % y = pdf(pd,x_values);
% % plot(x_values,y,'r','LineWidth',2)

end
%%
% burstMaxWidth = max([burstXwidth,burstYwidth],[],2) ;
% allVar = [burstDuration(1:end-1),burstArea(1:end-1),burstMaxWidth(1:end-1),...
%     peak(1:end-1),spaceDist,timeDist] ;
% [r,p] = corrcoef(allVar) ;
% % hist(areaBurst,20)
% % title(['Histogram of burst area, mean = ', meanArea(end)])
% % xlabel('Area (electrode^2)')
% 
% % print(gcf,[pwd,'/Results/random_noise/',dataFileName,...
% %        'HistArea_surr'],'-dpng')
% saveFileName = ['FullBurst',num2str(boundPrctile),...
%     '%',dataFileName,'_Start',int2str(timeStart),...
%     'minSize',num2str(minSize),'minTime',num2str(minBurstTime),'ms.mat'] ;
% save(saveFileName, 'burstDuration','burstXwidth','burstYwidth','burstArea'...
%     ,'peak','timeDist','spaceDist','r','p') ;    


%% make a movie
makeMovie = 1;
if makeMovie
    vidTitle = [pwd,'/Results/movies/New',dataFileName,num2str(timeStart2)] ;
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = 20 ;
    open(vidObj);
     fig=figure ;
     set(gcf,'Position',[260 23 1159 926])     
    
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

