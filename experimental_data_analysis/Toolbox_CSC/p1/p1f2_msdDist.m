function [superIdx,meanSlope] = p1f2_msdDist(WCentroids)
center = WCentroids ;
% % select bursts of duration more than 100ms
% center = [] ;
% count = 1 ;
% for iBurst = 1:length(Duration)
%     if Duration(iBurst)>100
%         center{count} = WCentroids{iBurst} ;
%         count = count+1 ;
%     end
% end

fullMSD = [] ;
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 1 ;
slopeAll = [] ;
superIdx = [] ;
% close all
for iBurst = 1:size(center,2)     % 198
    Trajectory = [] ;
   Trajectory = center{iBurst} ;
   Trajectory(:,3) = 1:size(center{iBurst},1) ;  
    
% for iCentre = 1:floor(size(center{iBurst} ,1)/4)
%     Trajectory(iCentre,:) = mean(center{iBurst}((iCentre-1)*4+1:4*iCentre,:)) ;
% end
% if length(1:4:size(center{iBurst},1)) > size(Trajectory,1)
%     Trajectory(:,3) = 1:4:size(center{iBurst},1)-4 ;
% else
%     Trajectory(:,3) = 1:4:size(center{iBurst},1) ;
% end
    itaQ = nan(5,1) ;
    countQ = 1 ;
    
    for powerQ = 2% 1:0.2:5
        % [MSD,tau] = get_MD(Trajectory,powerQ) ;%,grid_size);
        [MSD,tau] = get_MSD(Trajectory) ;%,grid_size);
    tempMSD = [MSD,tau] ;
    fullMSD = [fullMSD;tempMSD] ;
    
    p_temp = [];
    norErr = [] ;
    maxStart = 30 ;               % minimum tau used to fit MSD
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
    % discard the MSD with error rate more than 10%
    if (min(norErr)<0.1)
        bestIdx = find(norErr == min(norErr)) ;
        pBest = p_temp(bestIdx(end),:) ;
        y = exp(polyval(pBest,log(tau))) ;
        slopeAll(iBurst) = pBest(1) ;
    else
       continue
    end
    % keep the superdiffusive one to do average    
    %if pBest(1)>1
        msdSuper = [msdSuper;tempMSD] ;
        superIdx(countSuper) = iBurst ;
        countSuper = countSuper + 1; 
    %end
    %plot for single burst
%     loglog(tau,MSD,'.-','MarkerSize',12)
%     % plot(tau,MSD,'.-')
%     
%     hold on
%     loglog(tau,y)
%     % plot(tau,y)
%     title(['MSD within a burst (',num2str(size(center{iBurst},1)),' ms) '])%,...
%         %'with \alpha = ', num2str(pAll(1)),' tau = ',num2str(maxStart+bestIdx(end))])
%     xlabel('\tau (ms)')
%     ylabel('Mean Square Distance (mm)')
%     
%     str = {'p = ',num2str(pAll(1))};
%     text(max(tau)+5,max(y)+5,str)
%     xlim([1 10^3])
%     pause
%     close all
    % itaQ(countQ) = pAll(1) ;
    % countQ = countQ+1 ;
    end
    % plot(linspace(1,5,length(itaQ)),itaQ,'ro-')
    % pause
    % close all
end

% plot for average MSD and distribution
% figure;
% sigIn = msdSuper ;
% [MSDsort,~,MSDIdx] = unique(sigIn(:,2)) ;
% 
% for idx = 1:length(MSDsort)
%     meanMSD(idx) = mean(sigIn(find(MSDIdx==idx),1)) ;
%     stdMSD(idx) = std(sigIn(find(MSDIdx==idx),1)) ;
% end
% 
% % plot(MSDsort,meanMSD,'o')
% % errorbar(MSDsort,meanMSD,stdMSD)
% loglog(MSDsort,meanMSD)
% % set(gca,'yscale','log','xscale','log')
% hold on
% 
% tau = 1:200 ;
% p = polyfit(log(MSDsort(1:20)'),log(meanMSD(1:20)),1) ;
% y = exp(polyval(p,log(tau))) ;
% slope = p(1) ;
% loglog(tau,y)
% % title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
% xlabel('\tau (ms)')
% ylabel('Mean Square Distance (electrode)')
% 
% str = {'p = ',num2str(p(1))};
% text(max(tau)+5,max(y)+5,str)

% figure
slopeAll(slopeAll == 0) = [] ;
meanSlope = mean(slopeAll) ;
hist(slopeAll,20)
% title('distribution of MSD for my144')
xlabel('m.s.d.')
ylabel('count')
