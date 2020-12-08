
%% Theta
movieSig =  sigIn{8} ; % 13
    % sigMin = min(movieSig(:)) ;
    sigMin = 1*std(movieSig(:)) ;
    % sigMax = max(movieSig(:))*0.4 ;
    sigMax = 2*std(movieSig(:)) + mean(movieSig(:)) ;
load('/import/headnode2/xlon3884/Results_data/Project2/BurstDetection9/08FullBurst95%_ma027_032_Band5.6569HzDecember18_16:15_v4.mat')
countIdx = zeros(size(WCentroids)) ; 
%%
subplot(1,2,1)
for burstNum = 83%:length(WCentroids)   
 timeCount = 0 ;
    for iTime = rangeFrame(burstNum,1)-5: rangeFrame(burstNum,1)+160%timeStart+120:timeEnd
    % timeStart = 60400 ; 34560; 5140; 
    % for iTime = timeStart:60600
        imagesc(linspace(0,10,size(movieSig,1)),linspace(0,10,size(movieSig,1)),movieSig(:,:,iTime))
    %    colorbar
    xlabel('x')
    ylabel('y')
    %    title(['Gamma oscillation at ', num2str(iTime-timeStart,'%3.0f'),' ms with patterns (black) and center (red)'])
        caxis([sigMin sigMax])
         
%         hold on
%         pattIdx = find(fullBinary(:,:,iTime)==1) ;
%         [x,y] = ind2sub(size(squeeze(fullBinary(:,:,iTime))), pattIdx) ;
%         plot(y/size(sigIn,1)*10,x/size(sigIn,1)*10,'k.','markersize',2)
idx = find(iTime>rangeFrame(:,1)&iTime<rangeFrame(:,2)) ; 
for iBurst = 1:length(idx)
    [I,J,K] = ind2sub(size(movieSig),patternIdx{idx(iBurst)}) ;
    tempJ = J(K==iTime) ;
    tempI = I(K==iTime) ;
    tempJ2= tempJ ;
    tempJ = tempJ(intersect(find(mod(tempJ,size(movieSig,1)/10/2) ==0) , find(mod(tempI,size(movieSig,1)/10/2)==0) )) ;
    tempI = tempI(intersect(find(mod(tempJ2,size(movieSig,1)/10/2) ==0) , find(mod(tempI,size(movieSig,1)/10/2)==0) )) ;
    % plot(tempJ/size(movieSig,1)*10,tempI/size(movieSig,1)*10,'k.','markersize',12)
end

        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;

        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            if countIdx(iBurst) == 0
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(movieSig,1)*10,[centre(timeC,2)]/size(movieSig,1)*10,'r.-','markersize',8)
                countIdx(iBurst) = countIdx(iBurst) + 1 ;
            else
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(movieSig,1)*10,[centre(timeC,2)]/size(movieSig,1)*10,'r.-','markersize',8)
                %for iCentre = 1:countIdx(iBurst)
                    plot(centre(1:timeC,1)/size(movieSig,1)*10,centre(1:timeC,2)/size(movieSig,1)*10,...
                        'r.-','markersize',4)
                    
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
        % pause(0.1)
        timeCount = timeCount+1 ;
        if timeCount<36
        cla
        end
    % end
    end
    %pause
    %cla
%   close(vidObj);
end
 
%% Gamma
movieSig =  sigIn{13} ; % 13
    % sigMin = min(movieSig(:)) ;
    sigMin = 1*std(movieSig(:)) ;
    % sigMax = max(movieSig(:))*0.4 ;
    sigMax = 2*std(movieSig(:)) + mean(movieSig(:)) ;
load('/import/headnode2/xlon3884/Results_data/Project2/BurstDetection9/13FullBurst95%_ma027_032_Band32HzJanuary21_10:55_v4.mat')
countIdx = zeros(size(WCentroids)) ; 
%%
subplot(1,2,2)
for burstNum = 222%:length(WCentroids)   
 timeCount = 0 ;
    for iTime = rangeFrame(burstNum,1)-5: rangeFrame(burstNum,1)+45%timeStart+120:timeEnd
    % timeStart = 60400 ; 34560; 5140; 
    % for iTime = timeStart:60600
        imagesc(linspace(0,10,size(movieSig,1)),linspace(0,10,size(movieSig,1)),movieSig(:,:,iTime))
    %    colorbar
    xlabel('x')
    ylabel('y')
    %    title(['Gamma oscillation at ', num2str(iTime-timeStart,'%3.0f'),' ms with patterns (black) and center (red)'])
        caxis([sigMin sigMax])
         
%         hold on
%         pattIdx = find(fullBinary(:,:,iTime)==1) ;
%         [x,y] = ind2sub(size(squeeze(fullBinary(:,:,iTime))), pattIdx) ;
%         plot(y/size(sigIn,1)*10,x/size(sigIn,1)*10,'k.','markersize',2)
idx = find(iTime>rangeFrame(:,1)&iTime<rangeFrame(:,2)) ; 
for iBurst = 1:length(idx)
    [I,J,K] = ind2sub(size(movieSig),patternIdx{idx(iBurst)}) ;
    tempJ = J(K==iTime) ;
    tempI = I(K==iTime) ;
    tempJ2= tempJ ;
    tempJ = tempJ(intersect(find(mod(tempJ,size(movieSig,1)/10/2) ==0) , find(mod(tempI,size(movieSig,1)/10/2)==0) )) ;
    tempI = tempI(intersect(find(mod(tempJ2,size(movieSig,1)/10/2) ==0) , find(mod(tempI,size(movieSig,1)/10/2)==0) )) ;
    % plot(tempJ/size(movieSig,1)*10,tempI/size(movieSig,1)*10,'k.','markersize',12)
end

        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;

        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            if countIdx(iBurst) == 0
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(movieSig,1)*10,[centre(timeC,2)]/size(movieSig,1)*10,'r.-','markersize',8)
                countIdx(iBurst) = countIdx(iBurst) + 1 ;
            else
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(movieSig,1)*10,[centre(timeC,2)]/size(movieSig,1)*10,'r.-','markersize',8)
                %for iCentre = 1:countIdx(iBurst)
                    plot(centre(1:timeC,1)/size(movieSig,1)*10,centre(1:timeC,2)/size(movieSig,1)*10,...
                        'r.-','markersize',4)
                    
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
        % pause(0.1)
        timeCount = timeCount+1 ;
        if timeCount<36
        cla
        end
    % end
    end
    % pause
    % cla
%   close(vidObj);
end
 