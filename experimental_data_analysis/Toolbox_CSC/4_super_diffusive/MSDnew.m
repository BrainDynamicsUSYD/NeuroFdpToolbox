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
countSuper = 1 ;
slopeAll = [] ;
superIdx = [] ;
close all
for iBurst = 1:size(center,2)     % 24, 26, 216
    if Duration(iBurst)<100
        continue
    end
    Trajectory = [] ;
    Trajectory = center{iBurst}*4 ;
    Trajectory(:,3) = 1:size(center{iBurst},1) ;
    itaQ = nan(5,1) ;
    countQ = 1 ;
    
    for powerQ = 2% 1:0.2:5
        [MSD,tau] = get_MD(Trajectory,powerQ) ;%,grid_size);
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
    %if pBest(1)>1.1
        msdSuper = [msdSuper;tempMSD] ;
        superIdx(countSuper) = iBurst ;
        countSuper = countSuper + 1; 
    %end
    
    
%     if pBest(1)<1.2
%         continue
%     end
%     % plot for single burst
%     loglog(tau,MSD,'.-')
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
%     str = {'p = ',num2str(pBest(1))};
%     text(max(tau)+5,max(y)+5,str)
% 
%     pause
%     close all



%     itaQ(countQ) = pBest(1) ;
%     countQ = countQ+1 ;
    end
    % plot(linspace(1,5,length(itaQ)),itaQ,'ro-')
    % pause
    % close all
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
    % errorbar(MSDsort,meanMSD,stdMSD)
    loglog(MSDsort,meanMSD)
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
figure;
center = WCentroids ;
% deltaTime = [2,4,8,16,32] ; %,100]; %30 ;   
deltaTime = [64,32,16,8,4,2] ;
iTa = 0.45 ;
% iTa = 0.54 ;
% iTa = 0.45 ;
legendInfo = [] ;
count = 0 ;
for delta = 1:length(deltaTime)%:-1:1
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    DisplaceNorm = [] ;
for iBurst = 1:size(center,2)  % length(superIdx) % 
    % curBurst = superIdx(iBurst) ;
    curBurst = iBurst ;
    posCenter = center{curBurst} ;
    % DisplaceTemp = sqrt(sum((posCenter(deltaT+1 :end,:) - posCenter(1: end-deltaT,:)).^2,2)) ;
    DisplaceTemp = posCenter(deltaT+1 :end,1) - posCenter(1: end-deltaT,1) ;
    Displace = [Displace;DisplaceTemp] ;
    DisplaceTempNorm = DisplaceTemp/(deltaT^iTa) ;
    DisplaceNorm = [DisplaceNorm;DisplaceTempNorm] ;
end
count = count+1 ;
if count>1
    allDisplaceTemp = [DisplaceNorm; allDisplaceTemp] ;
else
    allDisplaceTemp = DisplaceNorm;
end
    
[n,x] = histcounts(DisplaceNorm,'normalization','pdf') ;
% loglog(x(2:end),n,'o')
semilogy(x(2:end),n,'o')
hold on
legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
end

fitVar = DisplaceNorm ; % allDisplaceTemp ;  
pd = fitdist(fitVar,'stable') ;
x = linspace(x(1),x(end)+3,1000) ;
y = pdf(pd,x) ;
hold on ; semilogy(x,y,'r-');

pd = fitdist(fitVar,'normal') ;
y = pdf(pd,x) ;
hold on ; 
% loglog(x,y,'k--');
semilogy(x,y,'k--')
legend('t = 30 ms','Gaussian')
legendInfo{end+1} = ['\alpha stable fit'] ;
legendInfo{end+1} = ['Gaussian fit'] ;
legend(legendInfo)
% xlabel('electrode')
xlabel(['displacement/\eta^{',num2str(iTa),'}'])
% ylabel(['displacement/\eta^{',num2str(iTa),'}'])
ylabel('probability')
% title('Displacement distribution')


%% fluctuation
close all
fullMSD = [] ;
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 1 ;
slopeAll = [] ;
center = WCentroids ;

for iBurst =  1:size(center,2)     % 24, 26, 216
    Trajectory = [] ;
    Trajectory = center{iBurst}*4 ;
    Trajectory(:,3) = 1:size(center{iBurst},1) ;
    itaQ = nan(5,1) ;
    countQ = 1 ;
    
    powerQ = 1 ;
    [MD,tau] = get_MD(Trajectory,powerQ) ;%,grid_size);
    powerQ = 2 ;
    [MSD,tau] = get_MD(Trajectory,powerQ) ;%,grid_size);
    
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
        superIdx(countSuper) = iBurst ;
        countSuper = countSuper + 1; 
    end
%     %plot for single burst
%     loglog(tau,MSD,'.-')
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
% 
%     % pause
%     close all
%     itaQ(countQ) = pAll(1) ;
%     
%     plot(linspace(1,5,length(itaQ)),itaQ,'ro-')
%     pause
%     close all
    fluc = sqrt(MSD - MD.^2) ;
    
    loglog(tau,fluc) 
    pause
    close all
end

