%% MSD calculation and displacement distribution
%
% This code can be used to calculate the MSD of the trajectory with an
% optimal maximal tau and discard the ones with more than 10% error rate.
%
% author: Xian Long,  supervisor: Pulin Gong
% June 21, 2018

%% select only 1 burst
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
countSuper = 0 ;
slopeAll = [] ;

for iBurst = 1:size(center,2)
    Trajectory = [] ;
    Trajectory = center{iBurst} ;
    Trajectory(:,3) = 1:size(center{iBurst},1) ;
    

    [MSD,tau] = get_MSD(Trajectory) ;%,grid_size);
    
    tempMSD = [MSD,tau] ;
    fullMSD = [fullMSD;tempMSD] ;
    
    p_temp = [];
    norErr = [] ;
    maxStart = 20 ;               % minimum tau used to fit MSD
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
    if pBest(1)>1.1
        msdSuper = [msdSuper;tempMSD] ;
        countSuper = countSuper + 1;
    end
%     %plot for single burst
%     loglog(tau,MSD,'.-')
%     % plot(tau,MSD,'.-')
%     
%     hold on
%     loglog(tau,y)
%     % plot(tau,y)
%     title(['MSD within a burst (',num2str(size(center{iBurst},1)),' ms) ',...
%         'with \alpha = ', num2str(pAll(1)),' tau = ',num2str(maxStart+bestIdx(end))])
%     xlabel('\tau (ms)')
%     ylabel('Mean Square Distance (electrode)')
%     
%     str = {'p = ',num2str(pAll(1))};
%     text(max(tau)+5,max(y)+5,str)
% 
%     pause
%     close all
end

%% plot for average MSD and distribution
figure;
sigIn = msdSuper ;
    [MSDsort,~,MSDIdx] = unique(sigIn(:,2)) ;
    
    for idx = 1:length(MSDsort)
        meanMSD(idx) = mean(sigIn(find(MSDIdx==idx),1)) ;
        stdMSD(idx) = std(sigIn(find(MSDIdx==idx),1)) ;
    end
    
    % plot(MSDsort,meanMSD,'o')
    errorbar(MSDsort,meanMSD,stdMSD)
    set(gca,'yscale','log','xscale','log')
    hold on
    
    p = polyfit(log(MSDsort(1:20)'),log(meanMSD(1:20)),1) ;
    y = exp(polyval(p,log(tau))) ;
    slope = p(1) ;
    loglog(tau,y)
    title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
    xlabel('\tau (ms)')
    ylabel('Mean Square Distance (electrode)')
    
    str = {'p = ',num2str(p(1))};
    text(max(tau)+5,max(y)+5,str)

figure
slopeAll(slopeAll == 0) = [] ;
hist(slopeAll,40)
title('distribution of MSD for my144')
xlabel('MSD')


%% displacement distribution
close all
center = WCentroids ;
deltaTime = 30 ;   
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
for iBurst = 1:size(center,2)
    posCenter = center{iBurst} ;
    DisplaceTemp = sum((posCenter(deltaT+1 :end,:) - posCenter(1: end-deltaT,:)).^2,2) ;
    Displace = [Displace;DisplaceTemp] ;
end
[n,x] = hist(Displace,240) ;
loglog(x,n/sum(n),'o')
hold on
legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
end
% legend(legendInfo)
% pd = fitdist(Displace,'stable')
pd = fitdist(Displace,'normal')
y = pdf(pd,x) ;
hold on ; loglog(x,y,'r--');
legend('t = 30 ms','Gaussian')

