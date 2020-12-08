function job_burst_pattern_Gaussian(arrayID)
cd ..
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
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;      
        
    case 2
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 0;
        subBand = [13 , 30] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 3
        dataFileName = 'my144_101' ;
        timeStart = 0 ;
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
        
    case 4
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;        % for movie
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 5
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [13 , 30] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 6
        dataFileName = 'my147_53' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        Fs = 1.0173e3 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
        
    case 7
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 8
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [13 , 30] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 9
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
        
    case 10
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 11
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [13 , 30] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 12
        dataFileName = 'ma025_03' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
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
%          randNumHalf =  2*pi*rand(size(sigOri,3)/2-1,1) ;
%          randNum = [0;randNumHalf;0;-flip(randNumHalf) ]' ;
%          for iChannel = setdiff(1:100,badChannels)
%              freqSig = fft(sigOriTemp(iChannel,:)) ;
%              absFreq = abs(freqSig) ;
%              phaseFreq = angle(freqSig) ;
%              reconSig = ifft(absFreq.*exp(1i*phaseFreq)) ;
%     
%              surSigTemp(iChannel,:) = ifft(absFreq.*exp(1i*(phaseFreq+randNum))) ;
%          end
      sigIn = randn(10,10,fix(300*fsTemporal)) ;
    
%    surSig = sigOri ;
    
    % find subband signals
%     bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,[],0) ;
%     % bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
%     hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
%     
%     if plotAmp
%         sigIn = abs(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:300*fsTemporal) ));
%     else
%         sigIn = angle(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:300*fsTemporal) ));
%     end
%     
%     clearvars bandpassSig hilbertSig
    

    %% Interpolation and smoothing
    
    close all
    plotLength = fix(300*fsTemporal) ;
    numChannel = 20 ;
    sigSmooth = zeros(numChannel,numChannel,plotLength) ;
    % timeStart = 5000 ;
    
    addpath(genpath([pwd,'/ToolNeuroPatt']))
    addpath(genpath([pwd,'/ToolOthers/nanconv']))
    
    resizeScale = 2 ;
    for iTime = 1:plotLength
        timeSlot = iTime+timeStart ;
        %Try smoothing
        filtWidth = 3;
        filtSigma = 0.6;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
        sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);
    end
    clearvars smoothTemp
% sigSmooth = sigIn ;
burstFlag = 1 ;
if burstFlag == 1
    % find 95 percentile
    if ~plotAmp
        sigBinary = zeros(size(sigSmooth)) ;
        sigBinary(sigSmooth>pi*5/6) = 1 ;
        sigPlot = sigBinary.*sigSmooth ;
    else
        boundPrctile = 95 ;
        prcBound = prctile(sigSmooth(:),boundPrctile) ;
        sigBinary = zeros(size(sigSmooth)) ;
        sigBinary(sigSmooth>prcBound) = 1;
        sigPlot = sigBinary.*sigSmooth ;
    end
else
    numChannels = size(sigSmooth,1)*size(sigSmooth,2) ;
    sigReshape = reshape(sigSmooth,numChannels,[]) ;
    GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,badChannels,...
    numChannels) ;
    sigBinary = zeros(size(sigReshape)) ;
    sigBinary(GammaBurstEvent.is_burst==1) = 1;
    sigBinary = reshape(sigBinary,size(sigSmooth,1),size(sigSmooth,2),[]) ;
    sigPlot = sigBinary.*sigSmooth ;
end
    
    
    %%
    % find burst in 3D
    CC = bwconncomp(sigBinary) ;                % 3D
    B1 = regionprops(CC,'BoundingBox');
    boundary1 = cat(1, B1.BoundingBox);
    Area = regionprops(CC,'Area') ;            % for computing total scale of bursts
    count = 1;                                 % for counting bursts
    Centroids = [] ;
    WCentroids = [] ;
    instantScale = [] ;
    instantPeakAmp = [] ;
    instantTotalPower = [] ;
    rangeFrame = [] ;
    width = [] ;
    % Amp = [] ;
    
    
    for iBurst = 1: size(CC.PixelIdxList,2)
        currentIdx = CC.PixelIdxList{iBurst} ;
        % Duration1(iBurst) = boundary(iBurst,6) ;    % should be similar to Duration2
        % Duration2: store all burst duration before threshold
        burstTimeEnd = floor((currentIdx(end)-1)/numChannel.^2) +1 ;
        burstTimeStart = floor((currentIdx(1)-1)/numChannel.^2) +1 ;
        Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;    
        
        % temporal threshold
        if Duration2(iBurst) < minBurstTime
            continue
        end
        % spatial threshold
        if boundary1(iBurst,4)<3 && boundary1(iBurst,5)<3
            continue
        end
        % Amp = [Amp; sigPlot(currentIdx)];
        
        % burst properties
        Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration        
        patternScale(count) = Area(iBurst).Area ;  % for total scale
        
        burstIdxTemp = zeros(size(sigPlot)) ;
        burstIdxTemp(currentIdx) = 1 ;
        currentBurst = sigPlot.*burstIdxTemp ;
        sumAmp(count) = sum(currentBurst(:)) ;     % sum of amplitude
        peakAmp(count) = max(currentBurst(:)) ;    % peak amplitude
        
        % loop through each time frame to study instaneous properties within burst
        timeCount = 1 ;
        for iTime = burstTimeStart:burstTimeEnd  
            instantPattern{count}(:,:,timeCount) = currentBurst(:,:,iTime) ;          
            instantBinary = burstIdxTemp(:,:,iTime) ;
            instantScale{count}(timeCount,:) = sum(instantBinary(:)) ;  % instant scale
            instantPeakAmp{count}(timeCount,:) = max(max(currentBurst(:,:,iTime))) ;
            instantTotalPower{count}(timeCount,:) = sum(sum(currentBurst(:,:,iTime).^2) ) ;
            S = regionprops(instantBinary,instantPattern{count}(:,:,timeCount),{'Centroid','WeightedCentroid'});
            % calculate 
            B = regionprops(instantBinary,'BoundingBox');
            boundary = cat(1, B.BoundingBox);
            width{count}(timeCount,:) = (boundary(:,3)+boundary(:,4))/2 ;
            
            Centroids{count}(timeCount,:) = cat(1, S.Centroid);
            WCentroids{count}(timeCount,:) = cat(1, S.WeightedCentroid) ;
            timeCount = timeCount + 1;
        end
        rangeFrame(count,:) = [burstTimeStart,burstTimeEnd] ; 
        % 
        if count>1
            firstCentroidsLoc = squeeze(WCentroids{count}(1,:)) ;
            % calculate the distance of centroids of two bursts
            distCent(count) = sqrt(sum((firstCentroidsLoc - lastCentroidsLoc).^2)) ;
            
            firstCentroidsTime = burstTimeStart ;
            % calculate the time interval between bursts
            centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
        end
        lastCentroidsLoc = squeeze(WCentroids{count}(end,:)) ;
        lastCentroidsTime = burstTimeEnd ;
        count = count+1 ;
    end

%%
% burstMaxWidth = max([burstXwidth,burstYwidth],[],2) ;
% allVar = [burstDuration(1:end-1),burstArea(1:end-1),burstMaxWidth(1:end-1),...
%     peak(1:end-1),spaceDist,timeDist] ;
% [r,p] = corrcoef(allVar) ;
% hist(areaBurst,20)
% title(['Histogram of burst area, mean = ', meanArea(end)])
% xlabel('Area (electrode^2)')

% print(gcf,[pwd,'/Results/random_noise/',dataFileName,...
%        'HistArea_surr'],'-dpng')

set(gca,'YDir','normal')
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;

saveData = 0;
if saveData
saveFileName = ['FullBurst90%',dataFileName,'_Band',...
   num2str(subBand(1)),'Hz',dateStr,'.mat'] ;
save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp','rangeFrame',...
    'patternScale','instantScale','Duration','Duration2','distCent','centInterval','instantPeakAmp','instantTotalPower','width') ;    
end
%%
plotMovie = 1 ;
if plotMovie
for iBurst = 1:count-1
    close all
    vidTitle = [pwd,'/Project1/Results/SurBurstMoviesSpaceSmoothLarge/',dataFileName,...
        num2str(minBurstTime),num2str(arrayID),num2str(iBurst,'%03d'),dateStr] ;
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = 20 ;
    open(vidObj);
    fig=figure ;
    set(gcf,'Position',[260 23 1159 926])
	timeCount = 0 ;
    
    for iTime = rangeFrame(iBurst,1):rangeFrame(iBurst,2)
	
        timeSlot = iTime+timeStart2 ;
        subplot(2,2,1)
        imagesc(sigSmooth(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Raw signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,2)
        timeCount = timeCount + 1 ;
        imagesc(sigPlot(:,:,iTime))
        set(gca,'YDir','normal')
        title(['95% signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,3)
        imagesc(instantPattern{iBurst}(:,:,timeCount))
        set(gca,'YDir','normal')
        title(['Detected burst at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,4)
        % imagesc(sigSmooth(:,:,iTime))
        currentPoint = WCentroids{iBurst}(timeCount,:) ;
        % nextPoint = largestCentroids(iTime+1,:) ;
        plot(currentPoint(:,1),currentPoint(:,2),'.','MarkerSize',20)
        xlim([0,numChannel])
        ylim([0,numChannel])
        title(['Detected burst Centre at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        colorbar

        % pause(0.01)
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        cla
    end
    close(vidObj);
    pause(0.01)

end
end

end


