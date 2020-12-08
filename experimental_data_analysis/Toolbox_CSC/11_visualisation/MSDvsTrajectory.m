fullMSD = [] ;
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 0 ;
slopeAll = [] ;

close all
for iBurst = 47%1:size(center,2)
    centre = WCentroids{iBurst} ;

    Trajectory = [] ;
    Trajectory = centre ;
    Trajectory(:,3) = 1:size(centre,1) ;
    

    [MSD,tau] = get_MSD(Trajectory) ;%,grid_size);
    
    tempMSD = [MSD,tau] ;
    fullMSD = [fullMSD;tempMSD] ;
    
    p_temp = [];
    norErr = [] ;
    maxStart = 20 ;
    tau_max = maxStart:min(100,size(MSD,1)) ;
    for fitIdx = 1:length(tau_max)       
        fitRange = (1:tau_max(fitIdx)) ;
        [pAll,S] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
        y = exp(polyval(pAll,log(tau(fitRange)))) ;
        errorRate = mean(abs((MSD(fitRange)-y(fitRange))./MSD(fitRange))) ;
        p_temp(fitIdx,:) = pAll ;
        % norErr(fitIdx) = S.normr ;       
        norErr(fitIdx) = errorRate ;
    end
    if (min(norErr)<0.1)
        bestIdx = find(norErr == min(norErr)) ;
        pBest = p_temp(bestIdx(end),:) ;
        y = exp(polyval(pBest,log(tau))) ;
        slopeAll(iBurst) = pBest(1) ;
    else
       continue
    end
        
    if pBest(1)>1.1
        msdSuper = [msdSuper;tempMSD] ;
        countSuper = countSuper + 1;
    end
    %plot
    loglog(tau,MSD,'.-')
    % plot(tau,MSD,'.-')
    
    hold on
    loglog(tau,y)
    % plot(tau,y)
    title(['MSD (',num2str(size(centre,1)),' ms) ',...
        'with \alpha = ', num2str(pAll(1)),' tau = ',num2str(maxStart+bestIdx(end))])
    xlabel('\tau (ms)')
    ylabel('Mean Square Distance (electrode)')
    
    str = {'p = ',num2str(pAll(1))};
    text(max(tau)+5,max(y)+5,str)

    
    figure;
    for j = 1: size(centre,2)
        %for i = 1:size(center,1)
        plot(centre(:,1),centre(:,2),'.-','markersize',20); hold on
        %end
        
        xlim([0,10])
        ylim([0,10])
        title('Trajectory of the weighted centre of the pattern')
    end
end
    
   




