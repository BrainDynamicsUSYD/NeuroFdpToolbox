function job_burst_pattern_freq(arrayID)
cd ..
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%

switch arrayID
    case 1
        dataFileName = 'ma025_03' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 80] ;         % Gamma, beta, ...
        subBand = [70, 80] ;          % Gamma and Delta
        subBand = [14, 21] ;          % all bands
        subBand = [1,3] ;             % all bands with log, equal BW
        subBand = [1,3] ;             % all bands with linear, equal BW
        subBand = [1,10] ;            % all bands linear, large BW
        subBand = [0.5,1.5] ;         % variable bands
        subBand = [0.75 1.25] ;       % log with different BW
        % subBand = [1 10] ;            % test
        
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;    
        
    case 2
        dataFileName = 'ma025_03' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [13 , 30] ;
        subBand = [50, 60] ;
        subBand = [11, 13] ;
        subBand = [3,5] ;
        subBand = [5,7] ;             % all bands with linear, equal BW
        subBand = [11,20] ;            % all bands linear, large BW
        subBand = [1.5,2.5] ;         % variable bands
        subBand = [1.5 2.5] ;         % log with different BW
        % subBand = [11 18] ;            % test
        
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 3
        dataFileName = 'ma025_03' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [8 , 13] ;
        subBand = [30, 40] ;
        subBand = [9, 10] ;
        subBand = [7,9] ;
        subBand = [9, 11] ;             % all bands with linear, equal BW
        subBand = [21,30] ;            % all bands linear, large BW
        subBand = [2.5,4] ;         % variable bands
        subBand = [3 5] ;         % log with different BW
        % subBand = [19 25] ;            % test
        
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 4
        dataFileName = 'ma025_03' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [4 , 8] ;
        subBand = [3, 4] ;
        subBand = [7, 8] ;
        subBand = [15,17] ;
        subBand = [13,15] ;             % all bands with linear, equal BW
        subBand = [31,40] ;            % all bands linear, large BW
        subBand = [4,6] ;         % variable bands
        subBand = [6 10] ;         % log with different BW
        % subBand = [26 30] ;            % test
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 5
        dataFileName = 'ma025_03' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        subBand = [1, 2] ;
        subBand = [5, 6] ;
        subBand = [31,33] ;
        subBand = [17,19] ;             % all bands with linear, equal BW
        subBand = [41,50] ;            % all bands linear, large BW
        subBand = [6,8] ;         % variable bands
        subBand = [12 20] ;         % log with different BW
        % subBand = [31 34] ;            % test
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 6
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [30 , 80] ;
        subBand = [70, 80] ;
        subBand = [63,65] ;
        subBand = [21,23] ;             % all bands with linear, equal BW
        subBand = [51,60] ;            % all bands linear, large BW
        subBand = [8,10] ;         % variable bands
        subBand = [24 40] ;         % log with different BW
        % subBand = [35 37] ;            % test
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 7
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [13 , 30] ;
        subBand = [50, 60] ;
        subBand = [127,129] ;
        subBand = [25,27] ;             % all bands with linear, equal BW
        subBand = [61,70] ;            % all bands linear, large BW
        subBand = [10,13] ;         % variable bands
        subBand = [48 80] ;         % log with different BW
        % subBand = [38 39] ;            % test
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 8
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [8 , 13] ;
        subBand = [30, 40] ;
        subBand = [29,31] ;             % all bands with linear, equal BW
        subBand = [71,80] ;            % all bands linear, large BW
        subBand = [13,20] ;         % variable bands
        subBand = [96 160] ;         % log with different BW
        % subBand = [40 41] ;            % test
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
        
    case 9
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [4 , 8] ;
        subBand = [3, 4] ;
        subBand = [33,35] ;             % all bands with linear, equal BW
        subBand = [81,90] ;            % all bands linear, large BW
        subBand = [20,30] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        
    case 10
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        subBand = [1, 2] ;
        subBand = [37,39] ;             % all bands with linear, equal BW
        subBand = [91,100] ;            % all bands linear, large BW
        subBand = [30,50] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
    case 11
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        subBand = [1, 2] ;
        subBand = [22, 29] ;
        subBand = [41,43] ;             % all bands with linear, equal BW
        subBand = [50,80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
        
     case 12
        dataFileName = 'ma027_032' ;
        timeStart = 0 ;     % for finding centres
        timeStart2 = 0;
        subBand = [0.5 , 4] ;
        subBand = [1, 2] ;
        subBand = [0.5, 80] ;
        subBand = [45,47] ;             % all bands with linear, equal BW
        subBand = [80,120] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
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
    
%    surSig = sigOri ;
    
    % find subband signals
    bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels,0) ;
    % bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
    hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
    
    if plotAmp
        sigIn = abs(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:300*fsTemporal) ));
    else
        sigIn = angle(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:300*fsTemporal) ));
    end
    
    clearvars bandpassSig hilbertSig
    

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
burstFlag = 0 ;
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
    GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,0,badChannels,...
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
        if boundary1(iBurst,4)<9 && boundary1(iBurst,5)<9
            continue
        end
        % Amp = [Amp; sigPlot(currentIdx)];
        
        % burst properties
        burstIdx(count) = iBurst ;
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
            S = regionprops(instantBinary,instantPattern{count}(:,:,timeCount),{'Centroid','WeightedCentroid'} );
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

% set(gca,'YDir','normal')
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;

saveData = 1;
if saveData
folderName = ['NeuroSpaceBurst'];
saveFileName = [folderName,'/FullBurstYifan',dataFileName,'_Band',...
   num2str(subBand(1)),'_',num2str(subBand(2)),'Hz',dateStr,'.mat'] ;
save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp','rangeFrame',...
    'patternScale','instantScale','Duration','Duration2','distCent','centInterval'...
    ,'instantPeakAmp','instantTotalPower','width','CC','burstIdx') ;    
end
%%
plotMovie = 0 ;
if plotMovie
for iBurst = 1:count-1
    close all
    vidTitle = [pwd,'/Project1/Results/burstMoviesSpaceSmoothLarge/',dataFileName,...
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
        timeCount = timeCount + 1 ;50
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

%%
plotBurst = 0 ;
if plotBurst
for iBurst = [17,23,39,46,49,51,52,64,768,769,780,788,789]
    close all
    if size(instantPattern{iBurst},3)<5*7+2
        continue
    end
    
    count = 1 ;
    figure;
    set(gcf,'Position',[680 401 1533 560]) ;
    for iTime = 2:5:5*7+2
        % iTime = rangeFrame(iBurst,1)+2:5:rangeFrame(iBurst,1)+5*7+2 % rangeFrame(iBurst,2)
        subplot(2,4,count)
        imagesc(instantPattern{iBurst}(:,:,iTime))
        set(gca,'YDir','normal')
%         title(['Detected burst at ',num2str(subBand),'Hz at ', ...
%             int2str(timeSlot),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
%         
%         imagesc(sigPlot(:,:,iTime))
%         set(gca,'YDir','normal')
%         % title(['95% signal at ',num2str(subBand),'Hz at ', ...
%         %    int2str(timeSlot/1000*1000),'ms'])
%         if plotAmp
%             caxis([0 0.5*max(sigSmooth(:))])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
        count = count+1 ;
        
    end
    savefig(['PatternPropagation',num2str(iBurst),'.fig'])
    
    %
    figure;
    set(gcf,'Position',[680 401 1533 560]) ;
    centre = WCentroids{iBurst} ;
    count = 1 ;
    for iTime = 2:1:5*7+2 % rangeFrame(iBurst,2)
        
        %subplot(2,4,count)
        plot([centre(iTime,1),centre(iTime+1,1)],[centre(iTime,2),...
            centre(iTime+1,2)],'k.-','markersize',20)

        count = count+1 ;
        xlim([0 20])
        ylim([0 20])
        hold on
    end
    savefig(['PatternTrajectory',num2str(iBurst),'.fig'])
    
end
end
% 