function pattDetection_v5(sigIn, sigBinary, params, flagMovie, flagSaveData,saveFileName)
% function for pattern detection
% in this version, 3D detection, max size of pattern are considered
%% spatial-temporal pattern detection based on continuity
CC = bwconncomp(sigBinary) ;               % continued points in space-time
B1 = regionprops(CC,'BoundingBox');        % bounding box of patterns
boundPatt = cat(1, B1.BoundingBox);
areaPatt = regionprops(CC,'Area') ;            % for computing total scale of patterns

%% initialization
countPatt = 1;                             % for counting patterns
Duration = [] ;                            % pattern duration
pattSize = [] ;                            % pattern size
Centroids = [] ;                           % geometry centre of patterns
WCentroids = [] ;                          % weighted centre of patterns
distCent = [] ;                            % distance between patterns
centInterval = [] ;                        % intervel between patterns
instantScale = [] ;                        % instantaneous size in patterns
instantPeakAmp = [] ;                      % instantaneous peak amplitudue in patterns
instantTotalPower = [] ;                   % instantaneous total power in patterns
rangeFrame = [] ;                          % the start and end time frames of patterns
width = [] ;                               % instantaneous size of patterns


% params.minpatternTime = 30 ;             % for 1 Gamma cycle (30Hz), 3 cycles (80)
sigPlot = sigBinary.*sigIn ;               % for calculating weighted centres
fullBinary = zeros(size(sigPlot)) ;        % store the final pattern index

%% further pattern detection
for iPatt = 1: size(CC.PixelIdxList,2)
    currentIdx = CC.PixelIdxList{iPatt} ;        % extract index for patterns
    DurationAll(iPatt) = boundPatt(iPatt,6) ; 

    pattTimeStart = boundPatt(iPatt,3)+0.5 ;     % note that the first frame
                                                 % starts with 0.5
    pattTimeEnd = pattTimeStart + boundPatt(iPatt,6) -1 ;

    
    % temporal threshold
    if DurationAll(iPatt) < params.minPattTime
        continue
    end
    % spatial threshold
    if boundPatt(iPatt,4)< params.minPattSize || boundPatt(iPatt,5)<params.minPattSize
        continue
    end
    % Amp = [Amp; sigPlot(currentIdx)];
    
    % pattern properties to be stored
    Duration(countPatt) = pattTimeEnd-pattTimeStart+1 ;    % duration
    pattSize(countPatt) = areaPatt(iPatt).Area ;  % for total scale
    
    pattIdxTemp = zeros(size(sigPlot)) ;
    pattIdxTemp(currentIdx) = 1 ;
    currentPatt = sigPlot.*pattIdxTemp ;
    % sumAmp(count) = sum(currentPatt(:)) ;     % sum of amplitude
    % peakAmp(count) = max(currentPatt(:)) ;    % peak amplitude
    
    % loop through each time frame to study instaneous properties within
    % patterns
    timeCount = 1 ;
    for iTime = pattTimeStart:pattTimeEnd
        % grab the current patterns
        instantBinary = pattIdxTemp(:,:,iTime) ;
        
        % find the max size pattern
        if timeCount == 1     % find the largest as patterns
            sizePatt = [] ;
            CC_temp = bwconncomp(instantBinary) ;                % 2D
            for i2DPatt = 1:CC_temp.NumObjects
                sizePatt(i2DPatt) = length(CC_temp.PixelIdxList{i2DPatt}) ;
            end
            [~,idxPatt] = max(sizePatt) ;
            instantBinary = zeros(size(instantBinary)) ;
            instantBinary(CC_temp.PixelIdxList{idxPatt}) = 1 ;
        else                   % find the next with the closest distance
            distTemp = [] ;
            centroid_temp = [] ;
            CC_temp = bwconncomp(instantBinary) ;                % 2D
            if length(CC_temp.PixelIdxList)>1
                S_temp = regionprops(CC_temp,'centroid') ;
                centroid_temp =  cat(1, S_temp.Centroid);
                lastCen = WCentroids{countPatt}(timeCount-1,:) ;
                for iPatt = 1:length(CC_temp.PixelIdxList)
                    distTemp(iPatt) = sum(bsxfun...
                        (@minus,centroid_temp(iPatt,:),lastCen).^2) ;
                end
                [~,minIdx] = min(distTemp) ;
                instantBinary = zeros(size(instantBinary)) ;
                instantBinary(CC_temp.PixelIdxList{minIdx}) = 1 ;
            end
        end
        fullBinary(:,:,iTime) = fullBinary(:,:,iTime) + instantBinary ;
        instantPattern{countPatt}(:,:,timeCount) = instantBinary.*sigPlot(:,:,iTime) ;
        
        instantScale{countPatt}(timeCount,:) = sum(instantBinary(:)) ;  % instant scale
        
        tempPeak = max(max(currentPatt(:,:,iTime))) ;
        instantPeakAmp{countPatt}(timeCount,:) = tempPeak ;
        
        instantTotalPower{countPatt}(timeCount,:) = sum(sum(currentPatt(:,:,iTime).^2) ) ;
        
        S = regionprops(instantBinary,instantPattern{countPatt}(:,:,timeCount),{'Centroid','WeightedCentroid'} );
        % calculate
        B2 = regionprops(instantBinary,'BoundingBox');
        
        boundPattInst = cat(1, B2.BoundingBox);
        width{countPatt}(timeCount,:) = (boundPattInst(3)+boundPattInst(4))/2 ;
        
        
        Centroids{countPatt}(timeCount,:) = cat(1, S.Centroid);
        WCentroids{countPatt}(timeCount,:) = cat(1, S.WeightedCentroid) ;
        timeCount = timeCount + 1;
    end
    rangeFrame(countPatt,:) = [pattTimeStart,pattTimeEnd] ;
    %
    if countPatt>1
        firstCentroidsLoc = squeeze(WCentroids{countPatt}(1,:)) ;
        % calculate the distance of centroids of two patterns
        distCent(countPatt) = sqrt(sum((firstCentroidsLoc - lastCentroidsLoc).^2)) ;
        
        firstCentroidsTime = pattTimeStart ;
        % calculate the time interval between patterns
        centInterval(countPatt) = firstCentroidsTime - lastCentroidsTime ;
    end
    lastCentroidsLoc = squeeze(WCentroids{countPatt}(end,:)) ;
    lastCentroidsTime = pattTimeEnd ;
    countPatt = countPatt+1 ;
end
CC_patterns = bwconncomp(fullBinary) ;       % continued points in space-time
patternIdx = CC_patterns.PixelIdxList ;      % save index to reduce size

%% saving data
if flagSaveData
    saveFileName = ['v5_',saveFileName] ;
    save(saveFileName, 'Centroids','WCentroids', 'rangeFrame',...
        'pattSize','instantScale','Duration','DurationAll','distCent',...
        'centInterval','instantPeakAmp','instantTotalPower','width','patternIdx') ;
end

%% sample movies
if flagMovie
% optional, load 'fullBinary' directly
% fullBinary = zeros(size(sigSmooth)) ;
% for iBurst = 1:length(patternIdx)
%     fullBinary(patternIdx{iBurst}) = 1;
% end
    movieSig = fullBinary.*sigIn ;
    sigMin = min(movieSig(:)) ;
    sigMax = max(movieSig(:))*0.8 ;
    timeStart = fix(5*1000)+1 ;
    timeEnd = fix(8*1000) ;
    
    % gif version
    frameCount = 0 ;
    im = [] ;
    fig = figure;
    set(gcf,'Position',[463 274 745 622]) ;
    for iTime = timeStart:timeEnd
        imagesc(movieSig(:,:,iTime))
        colorbar
        title(['time at', num2str(iTime-timeStart,'%3f')])
        caxis([sigMin sigMax])
        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;
        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            centre = WCentroids{iBurst} ;
            timeC = iTime- rangeFrame(iBurst,1)+1 ;
            plot([centre(timeC,1)],[centre(timeC,2)],'k.-','markersize',20)
        end
        xlim([0.5,size(sigIn,1)+0.5])
        ylim([0.5,size(sigIn,1)+0.5])
        
        % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        frameCount = frameCount+1 ;
        frame = getframe(fig);
        im{frameCount} = frame2im(frame);
        cla;
    end
    
    filename = [saveFileName,'.gif'];
    for idx = 1:frameCount
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
    %%
    % avi version
    % vidTitle = [pwd,'/experiment'] ;
    vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 20 ;
    open(vidObj);
    
    frameCount = 0 ;
    im = [] ;
    
    fig = figure;
    set(gcf,'Position',[78 279 1728 385]) ;
    % set(gcf,'Position',[375 49 1073 896]) ;  2,2

    for iTime = timeStart+60:timeEnd
        subplot(1,3,1)
        imagesc(sigIn(:,:,iTime))
        colorbar
        title(['original signals at ', num2str(iTime-timeStart,'%3f'),' ms'])
        caxis([sigMin sigMax])
        
        subplot(1,3,2)
        imagesc(sigBinary(:,:,iTime))
        colorbar
        title(['after threshold at ', num2str(iTime-timeStart,'%3f'),' ms'])
        caxis([0 1])
        
        subplot(1,3,3)
        imagesc(movieSig(:,:,iTime))
        colorbar
        title(['pattern detected at ', num2str(iTime-timeStart,'%3f'),' ms'])
        caxis([sigMin sigMax])
        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;
        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            centre = WCentroids{iBurst} ;
            timeC = iTime- rangeFrame(iBurst,1)+1 ;
            plot([centre(timeC,1)],[centre(timeC,2)],'k.-','markersize',20)
        end
        xlim([0.5,size(sigIn,1)+0.5])
        ylim([0.5,size(sigIn,1)+0.5])
        % subplot(2,2,4)

        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        cla
    end
    close(vidObj);

end
end