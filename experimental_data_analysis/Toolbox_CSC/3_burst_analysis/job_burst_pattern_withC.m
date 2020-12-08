function job_burst_pattern_withC(arrayID)
% for visualise origianl signals and burst, with centre
cd ..
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%
arrayID = 6 ;
switch arrayID
    case 1
        dataFileName = 'ma027_032'  ;
        timeStart = 600 ;     % for finding centres
        timeEnd  = 617 ;
        timeStart2 = 0;        % for movie
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;      
            resizeScale = 5 ;

    case 2
        dataFileName = 'ma027_032' ;
        timeStart = 2300 ;
        timeEnd = 2330;
        timeStart2 = 0 ;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
            resizeScale = 5 ;

    case 3
        dataFileName = 'ma027_032'  ;
        timeStart = 5140 ;
        timeEnd  = 5217 ;
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
        resizeScale = 5 ;
        
    case 4
        dataFileName = 'ma027_032' ;
        timeStart = 5366 ;     % for finding centres
        timeEnd  = 5444 ;
        timeStart2 = 0;        % for movie
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
            resizeScale = 5 ;

    case 5
        dataFileName = 'ma027_032'  ;
        timeStart = 5850 ;     % for finding centres
        timeEnd  = 5884 ;
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
            resizeScale = 5 ;

    case 6
        dataFileName = 'ma027_032'  ;
         timeStart = 2300 ;
        timeEnd = 2330;
                timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
            resizeScale = 10 ;

    case 7
        dataFileName = 'ma027_032' ;
        timeStart = 600 ;     % for finding centres
        timeEnd  = 617 ;
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
            resizeScale = 10 ;

    case 8
        dataFileName = 'ma027_032' ;
        timeStart = 5140 ;
        timeEnd  = 5217 ;
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4 ;
            resizeScale = 10 ;

    case 9
        dataFileName = 'ma027_032' ;
        timeStart = 5366 ;     % for finding centres
        timeEnd  = 5444 ;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 1 ;
            resizeScale = 10 ;

    case 10
        dataFileName = 'ma027_032'  ;
        timeStart = 5850 ;     % for finding centres
        timeEnd  = 5884 ;
        timeStart2 = 0;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        % totalSur = 200 ;
        minBurstTime = 30 ;
        % minSize = 4;
        resizeScale = 10 ;

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
%       surSig = reshape(real(surSigTemp),10,10,[]) ;
    
    surSig = sigOri ;
    
    % find subband signals
    bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels,0) ;
    % bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
    hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
    
    if plotAmp
        sigIn = abs(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:10*fsTemporal) ));
    else
        sigIn = angle(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:10*fsTemporal) ));
    end
    
    clearvars bandpassSig hilbertSig
    

    %% Interpolation and smoothing
    
    close all
    plotLength = fix(10*fsTemporal) ;
    numChannel = 10*resizeScale ;
    sigSmooth = zeros(numChannel,numChannel,plotLength) ;
    % timeStart = 5000 ;
    
    addpath(genpath([pwd,'/ToolNeuroPatt']))
    addpath(genpath([pwd,'/ToolOthers/nanconv']))
    
    for iTime = 1:plotLength
        timeSlot = iTime ;
        %Try smoothing
        filtWidth = 3;
        filtSigma = 0.6;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
        sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);
    end
    clearvars smoothTemp
    %%
% sigSmooth = sigIn ;
burstFlag = 0 ;
if burstFlag == 1
    % find 95 percentile
    if ~plotAmp
        sigBinary = zeros(size(sigSmooth)) ;
        sigBinary(sigSmooth>pi*5/6) = 1 ;
        sigPlot = sigBinary.*sigSmooth ;
    else
        boundPrctile = 99 ;
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
    instantPattern = [] ;
    
    for iBurst = 1: 1000%size(CC.PixelIdxList,2)
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

saveData = 0;
if saveData
saveFileName = ['FullBurst90%',dataFileName,'_Band',...
   num2str(subBand(1)),'Hz',dateStr,'.mat'] ;
save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp','rangeFrame',...
    'patternScale','instantScale','Duration','Duration2','distCent','centInterval','instantPeakAmp','instantTotalPower','width') ;    
end
%%
patternPart = zeros(numChannel,numChannel,270000) ;
    
for iBurst = 1:size(rangeFrame,1)
    patternPart(:,:,rangeFrame(iBurst,1):rangeFrame(iBurst,2)) = patternPart(:,:,...
    rangeFrame(iBurst,1):rangeFrame(iBurst,2))+ instantPattern{iBurst} ;
end
%%
plotMovie = 1 ;
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;
if plotMovie
    close all
%     vidTitle = [pwd,'/Results/Project2/Results/temp/',dataFileName,...
%         'Rescale_',num2str(resizeScale),'_',num2str(arrayID),dateStr] ;
%     vidObj = VideoWriter(vidTitle);
%     vidObj.FrameRate = 10 ;
%     open(vidObj);
    fig=figure ;
     set(gcf,'Position',[675 260 1910 701])
% 	timeCount = 0 ;
    oldWC = zeros(1,2) ;
    WC = zeros(1,2) ;
    num = 41 ;     % 47,64,65,25
    for iTime = rangeFrame(num,1)-5:rangeFrame(num,2)+5
        timeSlot = iTime+timeStart2 ;
        subplot(1,2,1)
        imagesc(sigSmooth(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Raw Gamma signals at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.3*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
%         subplot(2,2,2)
         
%         imagesc(sigPlot(:,:,iTime))
%         set(gca,'YDir','normal')
%         title(['95% signal at ',num2str(subBand),'Hz at ', ...
%             int2str(timeSlot/1000*1000),'ms'])
%         if plotAmp
%             caxis([0 0.5*max(sigSmooth(:))])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
        
        subplot(1,2,2)
%         axis equal;
%         hold on;
        h1 = imagesc(patternPart(:,:,iTime)) ;
        set(gca,'YDir','normal')
        title(['Detected burst at ', ...
            int2str(timeSlot),'ms (black dot: Amp. weighted centre)'])
        if plotAmp
            caxis([0 0.3*max(sigSmooth(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        uistack(h1,'bottom')
        
        hold on 
        tempImage = patternPart(:,:,iTime) ;
        tempImage(tempImage~=0) = 1 ;
        vBW = bwconncomp(tempImage) ;
        vIdx = getfield(vBW,'PixelIdxList') ;
        tempLength = [] ;
        if length(vIdx) ~= 0
            for i = 1:length(vIdx)
                tempLength(i) = length(vIdx{i}) ;
            end
            [~,idx] = max(tempLength) ;
            tempImage = zeros(numChannel,numChannel) ;
            tempImage(vIdx{idx}) = 1 ;
            S = regionprops(tempImage,patternPart(:,:,iTime),{'Centroid','WeightedCentroid'} );
            WC = getfield(S,'WeightedCentroid') ;
            h2 = plot(WC(1),WC(2),'k.','markersize',12) ;
        else
            if oldWC(1)~=0
                % delete(findobj(gca,'Type','line','Color','k'));
                oldWC = [0,0] ;
            end
        end
        if oldWC(1)~=0
            hold on
            h3 = plot([oldWC(1),WC(1)],[oldWC(2),WC(2)],'k') ;
        end
        oldWC = WC ;
%         subplot(2,2,4)
%         % imagesc(sigSmooth(:,:,iTime))
%         currentPoint = WCentroids{iBurst}(timeCount,:) ;
%         % nextPoint = largestCentroids(iTime+1,:) ;
%         plot(currentPoint(:,1),currentPoint(:,2),'.','MarkerSize',20)
%         xlim([0,numChannel])
%         ylim([0,numChannel])
%         title(['Detected burst Centre at ',num2str(subBand),'Hz at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         colorbar
% 
         pause(0.01)
 %       writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
         delete(h1)
         if exist('h2')
            delete(h2)
         end
    end
%     close(vidObj);
%     pause(0.01)

end


end


%%
loadBurstStat = 0 ;
if loadBurstStat
   load([pwd,'/Results_data/Project1/95%Yifan/FullBurstYifan2SDma027_032_Band30HzAugust14_12:41.mat']) 
end

%% for visualise specific burst
plotBurst = 0 ;
if plotBurst
for iBurst = 1
    % close all
    if size(instantPattern{iBurst},3)<5*7+2
        continue
    end
    
    figure;
    set(gcf,'Position',[680 401 1533 560]) ;
    for iTime = 1:100
        subplot(1,2,1)
        imagesc(sigIn)
        
        
        % iTime = rangeFrame(iBurst,1)+2:5:rangeFrame(iBurst,1)+5*7+2 % rangeFrame(iBurst,2)
        subplot(1,2,1)
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
    % savefig(['PatternPropagation',num2str(iBurst),'.fig'])
    
    %
    hold on 
    centre = WCentroids{iBurst} ;
    count = 1 ;
    for iTime = 2:1:5*7+2 % rangeFrame(iBurst,2)
        
        %subplot(2,4,count)
        plot([centre(iTime,1),centre(iTime+1,1)],[centre(iTime,2),...
            centre(iTime+1,2)],'k.-','markersize',20)

        count = count+1 ;
        xlim([0 numChannel])
        ylim([0 numChannel])
        hold on
    end
    % savefig(['PatternTrajectory',num2str(iBurst),'.fig'])
    
end
end
% 