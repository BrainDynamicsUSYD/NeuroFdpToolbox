sigIn = sigSmooth ;
    movieSig = sigSmooth ;
    sigMin = 0 ;
    sigMax = mean(sigSmooth(:)) + 3*std(sigSmooth(:)) ;
%     timeStart = fix(5*1000)+140 ;
%     timeEnd = fix(5*1000)+460 ;
    
%     fullBinary = zeros(size(sigSmooth)) ;
%     for iBurst = 1:length(patternIdx)
%     [I,J,K] = ind2sub(size(sigSmooth),patternIdx{iBurst}) ;
%     fullBinary(I,J,K) = 1 ;
%     end

timeStart = fix(32.2*1000)+180 ;    % fix(5*1000)+140 ;
    timeEnd = fix(64.2*1000)+460 ;  % fix(5*1000)+460 ;
fig = figure;
    % set(gcf,'Position',[78 279 1728 385]) ;  % 1,3
    % set(gcf,'Position',[375 49 1073 896]) ;  % 2,2
    % set(gcf,'Position',[638 225 687 605]) ;
    set(gcf,'Position',[514 184 828 700]) ;
 %   countIdx = zeros(size(WCentroids)) ;
    
    
    for iTime = 5700:15700+1000%timeStart+120:timeEnd
    % timeStart = 60400 ; 34560; 5140; 
    % for iTime = timeStart:60600
        imagesc(linspace(0,10,size(sigIn,1)),linspace(0,10,size(sigIn,1)),sigIn(:,:,iTime))
        colorbar
        title(['Gamma oscillation at ', num2str(iTime-timeStart,'%3.0f'),' ms with patterns (black) and center (red)'])
        caxis([sigMin sigMax])
         
%         hold on
%         pattIdx = find(fullBinary(:,:,iTime)==1) ;
%         [x,y] = ind2sub(size(squeeze(fullBinary(:,:,iTime))), pattIdx) ;
%         plot(y/size(sigIn,1)*10,x/size(sigIn,1)*10,'k.','markersize',2)
% idx = find(iTime>rangeFrame(:,1)&iTime<rangeFrame(:,2)) ; 
% for iBurst = 1:length(idx)
%     [I,J,K] = ind2sub(size(sigSmooth),patternIdx{idx(iBurst)}) ;
%     plot(J/size(sigSmooth,1)*10-0.5,I/size(sigSmooth,1)*10,'k.','markersize',2)
% end

%         hold on
%         burstIdx = find( iTime>rangeFrame(:,1) &iTime<rangeFrame(:,2)) ;
% 
%         for iBurst = burstIdx:burstIdx+length(burstIdx)-1
%             if countIdx(iBurst) == 0
%                 centre = WCentroids{iBurst} ;
%                 timeC = iTime- rangeFrame(iBurst,1)+1 ;
%                 plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',20)
%                 countIdx(iBurst) = countIdx(iBurst) + 1 ;
%             else
%                 centre = WCentroids{iBurst} ;
%                 timeC = iTime- rangeFrame(iBurst,1)+1 ;
%                 plot([centre(timeC,1)]/size(sigIn,1)*10,[centre(timeC,2)]/size(sigIn,1)*10,'r.-','markersize',20)
%                 %for iCentre = 1:countIdx(iBurst)
%                     plot(centre(1:timeC,1)/size(sigIn,1)*10,centre(1:timeC,2)/size(sigIn,1)*10,...
%                         'r--','markersize',8)
%                     
%                 %end
%                 countIdx(iBurst) = countIdx(iBurst) + 1 ;
%                 
%             end
%         end
% %         xlim([0.5,size(sigIn,1)+0.5])
% %         ylim([0.5,size(sigIn,1)+0.5])
%         xlim([0,10])
%         ylim([0,10])
        % subplot(2,2,4)

        % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        pause(0.1)
        cla
    end
%   close(vidObj);
   

