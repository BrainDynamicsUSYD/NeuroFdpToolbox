function job_burst_pattern_3D(arrayID)

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
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;      
        
    case 2
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 2;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 3
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 4;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 1 ;
        
    case 4
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;        % for movie
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;
        
    case 5
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 6
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 4;
        subBand = [30 , 80] ;
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
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4;
        
    case 8
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 9
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 4;
        subBand = [30 , 80] ;
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
        timeStart2 = 2;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        totalSur = 200 ;
        minBurstTime = 30 ;
        minSize = 4 ;
        
    case 12
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 4;
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
% for nSur = 1:totalSur
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


%% Interpolation and smoothing

close all
% plotLength = fix(100*fsTemporal) ;
sigSmooth = zeros(20,20,size(sigIn,3)) ;
% timeStart = 5000 ;

addpath(genpath([pwd,'/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))

resizeScale = 2 ;
for iTime = 1:size(sigIn,3)% plotLength
    timeSlot = iTime+timeStart ;
     %Try smoothing
    filtWidth = 3;
    filtSigma = 0.6;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);   
end
clearvars smoothTemp


%% find 95 percentile
sigSmoothTemp = sigIn ;
if ~plotAmp
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>pi*5/6) = 1 ;
    sigPlot = sigBinary.*sigSmoothTemp ;
else
    boundPrctile = 95 ;
    prcBound = prctile(sigSmoothTemp(:),boundPrctile) ;
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>prcBound) = 1;
    sigPlot = sigBinary.*sigSmoothTemp ;
end

%%
    CC = bwconncomp(sigBinary) ;                % 3D
    B = regionprops(CC,'BoundingBox');
    boundary = cat(1, B.BoundingBox);
    Area = regionprops(CC,'Area') ;
    count = 1;                                 % for counting bursts
    Centroids = [] ;
    WCentroids = [] ;
    instantScale = [] ;
    instantPeakAmp = [] ;
    instantTotalPower = [] ;
    rangeFrame = [] ;
    % Amp = [] ;
    
    burstIdxOld = zeros(size(sigPlot)) ;
    WCentroidsPlot = cell(size(sigBinary,3),1) ;
    %burstPlot = zeros(size(sigIn)) ;
    
    for iBurst = 1: size(CC.PixelIdxList,2)
        currentIdx = CC.PixelIdxList{iBurst} ;
        % Duration1(iBurst) = boundary(iBurst,6) ;    % should be similar to
        % Duration2
        burstTimeEnd = floor((currentIdx(end)-1)/100) +1 ;
        burstTimeStart = floor((currentIdx(1)-1)/100) +1 ;
        Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;    %
        if Duration2(iBurst) < minBurstTime
            continue
        end
        % Amp = [Amp; sigPlot(currentIdx)];
        
        Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration
        
        patternScale(count) = Area(iBurst).Area ;  % 4 total scale
        
        burstIdxTemp = zeros(size(sigPlot)) ;
        burstIdxTemp(currentIdx) = 1 ;
        currentBurst = sigPlot.*burstIdxTemp ;
        sumAmp(count) = sum(currentBurst(:)) ;     % sum of amplitude
        peakAmp(count) = max(currentBurst(:)) ;    % peak amplitude
      % 
        burstIdx = burstIdxTemp+burstIdxOld ;
        burstIdxOld = burstIdx;
        
        timeCount = 1 ;
        
        for iTime = burstTimeStart:burstTimeEnd  
            instantPattern = currentBurst(:,:,iTime) ; 
            % burstPlot(:,:,iTime) = burstPlot(:,:,iTime)+instantPattern ;
            instantBinary = burstIdxTemp(:,:,iTime) ;
            instantScale{count}(timeCount,:) = sum(instantBinary(:)) ;  % instant scale
            instantPeakAmp{count}(timeCount,:) = max(max(currentBurst(:,:,iTime))) ;
            instantTotalPower{count}(timeCount,:) = sum(sum(currentBurst(:,:,iTime).^2) ) ;
            S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});
            Centroids{count}(timeCount,:) = cat(1, S.Centroid);
            WCentroids{count}(timeCount,:) = cat(1, S.WeightedCentroid) ;
            timeCount = timeCount + 1;
            
            WCentroidsPlot{iTime} = [WCentroidsPlot{iTime},cat(1, S.WeightedCentroid)] ;
        end
        rangeFrame(count,:) = [burstTimeStart,burstTimeEnd] ; 
        if count>1
            firstCentroidsLoc = squeeze(WCentroids{count}(1,:)) ;
            distCent(count) = sqrt(sum((firstCentroidsLoc - lastCentroidsLoc).^2)) ;
            
            firstCentroidsTime = burstTimeStart ;
            centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
        end
        lastCentroidsLoc = squeeze(WCentroids{count}(end,:)) ;
        lastCentroidsTime = burstTimeEnd ;
        count = count+1 ;
        
    end


%% make a movie
sigBurst = sigIn.*burstIdx ;
makeMovie = 1;
if makeMovie
    vidTitle = [pwd,'/Results/movies/New',dataFileName,'_', num2str(boundPrctile) ,'%_',...
        num2str(minBurstTime),'ms',num2str(timeStart2),'.avi'] ;
    vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 20 ;
    open(vidObj);
     fig=figure ;
     set(gcf,'Position',[260 23 1159 926])     
    
    for iTime = 1:2000
        timeSlot = iTime+timeStart2 ;
        subplot(2,2,1)
        imagesc(sigIn(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        title(['Raw signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        if plotAmp
            caxis([0 0.6*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,4)
        % imagesc(sigSmooth(:,:,iTime+timeStart))       
        imagesc(sigBurst(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        title(['Detected burst at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        if plotAmp
            caxis([0 0.6*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,3)
        imagesc(sigPlot(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        title(['95% threshold at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        if plotAmp
            caxis([0 0.6*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,2)
        % imagesc(sigSmooth(:,:,iTime))
        currentPoint = WCentroidsPlot{iTime+timeStart} ;
        currentPointX = currentPoint(1:2:end) ;
        currentPointY = currentPoint(2:2:end) ;
        % nextPoint = largestCentroids(iTime+1,:) ;
        plot(currentPointX,currentPointY,'.','MarkerSize',20)
        xlim([0.5,10.5])
        ylim([0.5,10.5])
        title(['Detected burst Centre at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        colorbar
        currentPointX = [] ;
        currentPointY = [] ;
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        % pause(0.01)
        cla
        
    end
    close(vidObj);
end

