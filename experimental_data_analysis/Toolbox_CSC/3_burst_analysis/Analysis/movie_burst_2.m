function movie_burst(sigSmooth,WCentroids,patternIdx,rangeFrame,flagType)
sigIn = sigSmooth ;
    movieSig = sigSmooth ;
    % sigMin = min(movieSig(:)) ;
    sigMin = 0 ;
    % sigMax = max(movieSig(:))*0.4 ;
    sigMax = 3*std(sigSmooth(:)) + mean(sigSmooth(:)) ;
%     timeStart = fix(5*1000)+140 ;
%     timeEnd = fix(5*1000)+460 ;
    
%     fullBinary = zeros(size(sigSmooth)) ;
%     for iBurst = 1:length(patternIdx)
%     [I,J,K] = ind2sub(size(sigSmooth),patternIdx{iBurst}) ;
%     fullBinary(I,J,K) = 1 ;
%     end

timeStart = fix(32.2*1000)+180 ;    % fix(5*1000)+140 ;
    timeEnd = fix(64.2*1000)+460 ;  % fix(5*1000)+460 ;
    %%
    % close all

    if flagType == 1
        burstNum = 17 ;
    elseif flagType == 2
        burstNum = 409 ;
    end
    %%
    countIdx = zeros(size(WCentroids)) ;
    % close all
    % for burstNum = 420:440;  % 17,224, 43, 46, 97,98 ; sticky: 27,92,262,270,278,299,331,409
    % fig = figure;
    % set(gcf,'Position',[78 279 1728 385]) ;  % 1,3
    % set(gcf,'Position',[375 49 1073 896]) ;  % 2,2
    % set(gcf,'Position',[638 225 687 605]) ;
    % set(gcf,'Position',[514 184 828 700]) ;
    
    timeCount = 0 ;
    for iTime = rangeFrame(burstNum,1)-5: rangeFrame(burstNum,1)+30%timeStart+120:timeEnd
    % timeStart = 60400 ; 34560; 5140; 
    % for iTime = timeStart:60600
        imagesc(linspace(0,10,size(sigIn,1)),linspace(0,10,size(sigIn,1)),sigIn(:,:,iTime))
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
    [I,J,K] = ind2sub(size(sigSmooth),patternIdx{idx(iBurst)}) ;
    tempJ = J(K==iTime) ;
    tempI = I(K==iTime) ;
    tempJ2= tempJ ;
    tempJ = tempJ(intersect(find(mod(tempJ,size(sigSmooth,1)/10/2) ==0) , find(mod(tempI,size(sigSmooth,1)/10/2)==0) )) ;
    tempI = tempI(intersect(find(mod(tempJ2,size(sigSmooth,1)/10/2) ==0) , find(mod(tempI,size(sigSmooth,1)/10/2)==0) )) ;
    plot(tempJ/size(sigSmooth,1)*10,tempI/size(sigSmooth,1)*10,'k.','markersize',12)
end

        hold on
        burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;

        for iBurst = burstIdx:burstIdx+length(burstIdx)-1
            if countIdx(iBurst) == 0
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',8)
                countIdx(iBurst) = countIdx(iBurst) + 1 ;
            else
                centre = WCentroids{iBurst} ;
                timeC = iTime- rangeFrame(iBurst,1)+1 ;
                plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',8)
                %for iCentre = 1:countIdx(iBurst)
                    plot(centre(1:timeC,1)/size(sigIn,1)*10,centre(1:timeC,2)/size(sigIn,1)*10,...
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
%   close(vidObj);
   

