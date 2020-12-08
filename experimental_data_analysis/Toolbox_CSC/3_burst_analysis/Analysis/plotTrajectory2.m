for iBurst = 18%:size(patternIdx,2)
    close all
    temp = zeros(size(sigSmooth)) ;
    temp(patternIdx{iBurst}) = 1 ;
    temp2 = temp.*sigSmooth ;
    maxSig = max(sigSmooth(:)) ;
    minSig = 0 ;
    count = 1;
    timeStart = floor(patternIdx{iBurst}(1)/400)+10 ;
    for iTime =  timeStart:10:ceil(patternIdx{iBurst}(end)/400)
        subplot(1,4,count)
        imagesc(temp2(:,:,iTime))
        caxis([0 0.4*maxSig])
        set(gca,'YDir','normal');
        count = count+1;
    end
    disp(['Burst number = ',num2str(iBurst)])
    
    
    figure;
    center = WCentroids{iBurst} ;
    subplot(1,4,1)
   
    plot(center(1,1),center(1,2),'.-','markersize',20); hold on
    xlim([0,20])
    ylim([0,20])
    for iImage = 2 :4
        
        subplot(1,4,iImage)        
        %for iTime = 1+10*(iImage-2):10+10*(iImage-2)
        iTime = 1+10*(iImage-2);
        plot(center(iTime:iTime+10,1),center(iTime:iTime+10,2),'.-','markersize',20); hold on
        %end
        %end
        xlim([0,20])
        ylim([0,20])
        % title(num2str(j))

    end
end