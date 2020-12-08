function pattDetection_v7(sigIn, sigBinary, params, flagMovie, flagSaveData,saveFileName)
% function for pattern detection
% for project2
%% spatial-temporal pattern detection based on continuity
s = size(sigBinary);
    CC = bwconncomp(sigBinary,6) ;                % 3D
    % per = [1 1 0];
    % CC = CC2periodic(CC,per);1
    B1 = regionprops(CC,'BoundingBox');
    boundary1 = cat(1, B1.BoundingBox);
    
    Area = regionprops(CC,'Area') ;            % for computing total scale of bursts
    count = 1;                                 % for counting bursts
    numChannel = size(sigIn{iBand},1) ;
    minBurstTime = fsTemporal/cfreqso(iBand) ;         % 1 cycle in ms
    Centroids = [] ;
    WCentroids = [] ;
    instantScale = [] ;
    instantPeakAmp = [] ;
    instantTotalPower = [] ;
    rangeFrame = [] ;
    width = [] ;
    burstDetect = [] ;
    
    % Amp = [] ;
    num_Burst = 0 ;
    num_Burst_longtime = 0 ;
    w = size(sigBinary,2) ;
    close all
    
    for iBurst = 1: size(CC.PixelIdxList,2)
        currentIdx = CC.PixelIdxList{iBurst} ;
        % Duration1(iBurst) = boundary(iBurst,6) ;    % should be similar to Duration2
        % Duration2: store all burst duration before threshold
        burstTimeEnd = floor((currentIdx(end)-1)/numChannel.^2) +1 ;
        burstTimeStart = floor((currentIdx(1)-1)/numChannel.^2) +1 ;
        Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;
        num_Burst = num_Burst+1 ;
        
        if isinteger(iBurst*100/size(CC.PixelIdxList,1))
            fprintf(['burst at ',num2str(iBurst*100/size(CC.PixelIdxList,1)),'%'])
        end
        
        % temporal threshold
        if Duration2(iBurst) < minBurstTime
            % timelimit = 1
            continue
        end
        
        num_Burst_longtime = num_Burst_longtime+1 ;
        % spatial threshold
        if boundary1(iBurst,4)<4 && boundary1(iBurst,5)<4
            % spacelimit = 1
            continue
        end
        % Amp = [Amp; sigPlot(currentIdx)];
        
        % burst properties
        Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration
        patternScale(count) = Area(iBurst).Area ;  % for total scale
        
        burstIdxTemp = zeros(size(sigPlot{iBand})) ;
        burstIdxTemp(currentIdx) = 1 ;
        currentBurst = sigPlot{iBand}.*burstIdxTemp ;
        sumAmp(count) = sum(currentBurst(:)) ;     % sum of amplitude
        peakAmp(count) = max(currentBurst(:)) ;    % peak amplitude
        
        burstStart(count) = burstTimeStart ;
        burstDetect.PixelIdxList{count} = currentIdx ;
        
        %% loop through instantaneous frames to find instantaneous size
        timeCount = 1 ;
        for iTime = burstTimeStart:burstTimeEnd
            instantPattern{count}(:,:,timeCount) = currentBurst(:,:,iTime) ;
            
            % CC_2D = bwconncomp(burstIdxTemp(:,:,iTime)) ;    % 2D check how many patterns
            % num_pattern{count}(timeCount,:) = length(CC_2D.PixelIdxList) ;
            
            % some statistics
            instantBinary = burstIdxTemp(:,:,iTime) ;
            %         instantBinary_full = [instantBinary(w/2+1:end,w/2+1:end),instantBinary(w/2+1:end,1:w/2+1),...
            %             instantBinary(w/2+1:end,w/2+1:end),instantBinary(w/2+1:end,1:w/2+1);...
            %             instantBinary(1:w/2,w/2+1:end),instantBinary(1:w/2,1:w/2),...
            %             instantBinary(1:w/2,w/2+1:end),instantBinary(1:w/2,1:w/2);...
            %             instantBinary(w/2+1:end,w/2+1:end),instantBinary(w/2+1:end,1:w/2+1),...
            %             instantBinary(w/2+1:end,w/2+1:end),instantBinary(w/2+1:end,1:w/2+1);...
            %             instantBinary(1:w/2,w/2+1:end),instantBinary(1:w/2,1:w/2),...
            %             instantBinary(1:w/2,w/2+1:end),instantBinary(1:w/2,1:w/2)];
            
            CC_2D = bwconncomp(instantBinary) ;
            % per = [1 1];
            % CC_2D = CC2periodic(CC_2D,per);
            B1_2D = regionprops(CC_2D,'BoundingBox');
            S_2D = regionprops(CC_2D,'Centroid');
            Area_2D = regionprops(CC_2D,'Area') ;
            area_2D = cat(1,Area_2D.Area) ;
            num_pattern{count}(timeCount,:) = length(CC_2D.PixelIdxList) ;
            instantScale{count}(timeCount,:) = sum(area_2D) ;
            
            S = regionprops(instantBinary,instantPattern{countPatt}(:,:,timeCount),{'Centroid','WeightedCentroid'} );
            Centroids{count}(timeCount,:) = cat(1, S.Centroid);
            WCentroids{count}(timeCount,:) = cat(1, S.WeightedCentroid) ;
                    
            timeCount = timeCount + 1;

        end
        count = count+1 ;
    end
%% saving data
   end
    t = datetime('now') ;
    dateStr = datestr(t,'mmmmdd_HH:MM') ;
    
    if ~exist('sumAmp')
        disp(['no burst for ',num2str(iBand)])
        continue
    end
    
    disp(['saving ',num2str(iBand)])
    
    saveFileName1 = [num2str(iBand,'%02d'),'FullBurst95%',...
            num2str(cfreqso(iBand)),'Hz_',dataFileName,dateStr,num2str(arrayID+8,'%02d'),'.mat']
    
    saveData = 1;
    if saveData
        save([saveFolderName,saveFileName1],  'sumAmp','peakAmp','burstStart',...
            'patternScale','Duration','Duration2','num_Burst','burstDetect','CC','instantScale','WCentroids','Centroids') ;
    end

%% sample movies
if flagMovie
% optional, load 'fullBinary' directly
% fullBinary = zeros(size(sigSmooth)) ;
% for iBurst = 1:length(patternIdx)
%     fullBinary(patternIdx{iBurst}) = 1;
% end
    % movieSig = fullBinary.*sigIn ;
    movieSig = sigSmooth ;
    sigMin = min(movieSig(:)) ;
    sigMax = max(movieSig(:))*0.4 ;
    timeStart = fix(5*1000)+140 ;
    timeEnd = fix(5*1000)+460 ;
    
%     % gif version
%     frameCount = 0 ;
%     im = [] ;
%     fig = figure;
%     set(gcf,'Position',[463 274 745 622]) ;
%     for iTime = timeStart:timeEnd
%         imagesc(movieSig(:,:,iTime))
%         colorbar
%         title(['time at', num2str(iTime-timeStart,'%3f')])
%         caxis([sigMin sigMax])
%         hold on
%         burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;
%         for iBurst = burstIdx:burstIdx+length(burstIdx)-1
%             centre = WCentroids{iBurst} ;
%             timeC = iTime- rangeFrame(iBurst,1)+1 ;
%             plot([centre(timeC,1)],[centre(timeC,2)],'k.-','markersize',20)
%         end
%         xlim([0.5,size(sigIn,1)+0.5])
%         ylim([0.5,size(sigIn,1)+0.5])
%         
%         % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
%         frameCount = frameCount+1 ;
%         frame = getframe(fig);
%         im{frameCount} = frame2im(frame);
%         cla;
%     end
%     
%     filename = [saveFileName,'.gif'];
%     for idx = 1:frameCount
%         [A,map] = rgb2ind(im{idx},256);
%         if idx == 1
%             imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
%         else
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
%         end
%     end
    %%
%     % avi version
%     % vidTitle = [pwd,'/experiment'] ;

%     vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);
     
%     frameCount = 0 ;
%     im = [] ;
%     
%     fig = figure;
%     set(gcf,'Position',[78 279 1728 385]) ;  % 1,3
%     % set(gcf,'Position',[375 49 1073 896]) ;  % 2,2
%     % set(gcf,'Position',[638 225 687 605]) ;
% 
%     for iTime = timeStart+60:timeEnd
%         subplot(1,3,1)
%         imagesc(sigIn(:,:,iTime))
%         colorbar
%         title(['original signals at ', num2str(iTime-timeStart,'%3f'),' ms'])
%         caxis([sigMin sigMax])
%         
%         subplot(1,3,2)
%         imagesc(sigBinary(:,:,iTime))
%         colorbar
%         title(['after threshold at ', num2str(iTime-timeStart,'%3f'),' ms'])
%         caxis([0 1])
%         
%         subplot(1,3,3)
%         imagesc(movieSig(:,:,iTime))
%         colorbar
%         title(['pattern detected at ', num2str(iTime-timeStart,'%3f'),' ms'])
%         caxis([sigMin sigMax])
%         hold on
%         burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;
%         for iBurst = burstIdx:burstIdx+length(burstIdx)-1
%             centre = WCentroids{iBurst} ;
%             timeC = iTime- rangeFrame(iBurst,1)+1 ;
%             plot([centre(timeC,1)],[centre(timeC,2)],'k.-','markersize',20)
%         end
%         xlim([0.5,size(sigIn,1)+0.5])
%         ylim([0.5,size(sigIn,1)+0.5])
%         % subplot(2,2,4)
% 
%         writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
%         cla
%     end
%     close(vidObj);

    
    % avi version
    % vidTitle = [pwd,'/experiment'] ;
%     vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);
    
    frameCount = 0 ;
    im = [] ;
    
    fig = figure;
    % set(gcf,'Position',[78 279 1728 385]) ;  % 1,3
    % set(gcf,'Position',[375 49 1073 896]) ;  % 2,2
    % set(gcf,'Position',[638 225 687 605]) ;
    set(gcf,'Position',[514 184 828 700]) ;
    countIdx = zeros(size(WCentroids)) ;
    for iTime = timeStart+120:timeEnd
    % timeStart = 60400 ;
    % for iTime = timeStart:60600
        imagesc(linspace(0,10,size(sigIn,1)),linspace(0,10,size(sigIn,1)),sigIn(:,:,iTime))
        colorbar
        title(['Gamma oscillation at ', num2str(iTime-timeStart,'%3.0f'),' ms with patterns (black) and center (red)'])
        caxis([sigMin sigMax])
         
        hold on
        pattIdx = find(fullBinary(:,:,iTime)==1) ;
        [x,y] = ind2sub(size(squeeze(fullBinary(:,:,iTime))), pattIdx) ;
        plot(y/size(sigIn,1)*10,x/size(sigIn,1)*10,'k.','markersize',2)
        
        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;

        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            if countIdx(iBurst) == 0
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',20)
                countIdx(iBurst) = countIdx(iBurst) + 1 ;
            else
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',20)
                %for iCentre = 1:countIdx(iBurst)
                    plot(centre(1:timeC,1)/size(sigIn,1)*10,centre(1:timeC,2)/size(sigIn,1)*10,...
                        'r--','markersize',8)
                    
                %end
                countIdx(iBurst) = countIdx(iBurst) + 1 ;
                
            end
        end
%         xlim([0.5,size(sigIn,1)+0.5])
%         ylim([0.5,size(sigIn,1)+0.5])
        xlim([0,10])
        ylim([0,10])
        % subplot(2,2,4)

        % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        pause
        cla
    end
%   close(vidObj);
   

end
end