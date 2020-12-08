addpath(genpath([pwd,'/Toolbox_CSC']))
%% check how many bursts at a same time

burst = zeros(size(Centroids,1),1) ;

for iBurst = 1:size(Centroids,2)
    temp = zeros(size(Centroids,1),1) ;
    temp(Centroids{iBurst}(:,1)~=0)  = 1 ;
    burst = burst + temp ;
    
end
plot(burst)
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

for iBurst = 39%1:size(center,2)
    Trajectory = [] ;
    Trajectory = center{iBurst} ;
    Trajectory(:,3) = 1:size(center{iBurst},1) ;
    

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
%     if (min(norErr)<0.1)
%         bestIdx = find(norErr == min(norErr)) ;
%         pBest = p_temp(bestIdx(end),:) ;
%         y = exp(polyval(pBest,log(tau))) ;
%         slopeAll(iBurst) = pBest(1) ;
%     else
%        continue
%     end
%         
%     if pBest(1)>1.1
%         msdSuper = [msdSuper;tempMSD] ;
%         countSuper = countSuper + 1;
%     end
    %plot
    loglog(tau,MSD,'.-')
    % plot(tau,MSD,'.-')
    
    hold on
    loglog(tau,y)
    % plot(tau,y)
    title(['MSD within a burst (',num2str(size(center{iBurst},1)),' ms) ',...
        'with \alpha = ', num2str(pAll(1)),' tau = ',num2str(maxStart+bestIdx(end))])
    xlabel('\tau (ms)')
    ylabel('Mean Square Distance (electrode)')
    
    str = {'p = ',num2str(pAll(1))};
    text(max(tau)+5,max(y)+5,str)

    pause
    close all
end
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


%% distribution
msdTemp = msdSuper;
% msdTemp = msdSuper ;
[MSDsort,~,MSDIdx] = unique(msdTemp(:,2)) ;
close all;
tauSelect = [2,4,8,16,32,64] ;
eta = tauSelect.^1.21 ;
%eta = tauSelect.^1.5 ;
rouNormFull = [] ;
for idx = 1:6 
    % figure;
    rouTemp = msdTemp( find(MSDIdx==tauSelect(idx) ))/eta(idx) ;
    [rouHist,n] = hist(rouTemp, 400) ;
    rouNormFull = [rouNormFull;rouTemp] ;
    loglog(n,rouHist/sum(rouHist),'o')
    hold on
end
legend('1','2','3','4','5','6')
figure;
for idx = 1:6 
    % figure;
    [r,n] = hist(msdTemp(find(MSDIdx==tauSelect(idx))), 400) ;
    loglog(n,r/sum(r),'o')
    hold on
end
legend('1','2','3','4','5','6')

%%
[alpha, xmin, L] = plfit(rouNormFull) 
plplot(x, xmin, alpha, 'bo')

%%
close all
figure
xmin = 0.05 ;
sigIn = rouNormFull ;
sigIn = sigIn(sigIn>xmin) ;

pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
[lambdaHat1,lambdaCI] = mle(sigIn, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
L1 = sum(log(pf_power(sigIn,lambdaHat1 ))) ;

z = linspace(0.015,0.5,200);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_power(x,lambdaHat1))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat1));

% MLE of power law with cut off
figure;
pf_powerC = @(x,alpha,lambda) x.^-alpha.*exp(-lambda*x)*lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
[lambdaHat2,lambdaCI] = mle(sigIn, 'pdf',pf_powerC, 'start',[1,0], 'lowerbound',[1,0], 'upperbound', [3,10]) ;
L2 = sum(log(pf_powerC(sigIn,lambdaHat2(1), lambdaHat2(2) ))) ;


histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_powerC(x,lambdaHat2(1),lambdaHat2(2)))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_powerC(x, lambdaHat2(1),lambdaHat2(2)));


% MLE of power law with cut off
figure;
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 50) ;
L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;

histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_exp(x,lambdaHat3))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));

% %% get MD
% slope = [] ;
% fullMD = [] ;
% meanMD = [] ;
% stdMD = [] ;
% for iBurst = 1:size(Centroids,2)
%     Trajectory = [] ;
%     Trajectory = center{iBurst} ;
%     Trajectory(:,3) = 1:size(center{iBurst},1) ;
%     
% 
%     [MD,tau] = get_MD(Trajectory) ;%,grid_size);
%     
%     tempMD = [MD,tau] ;
%     fullMD = [fullMD;tempMD] ;
% 
%     
%     pAll = polyfit(log(tau(1:20)),log(MD(1:20)),1) ;
%     y = exp(polyval(pAll,log(tau))) ;
%     slopeAll(iBurst) = pAll(1) ;
%     
%     %plot
% %     loglog(tau,MSD,'.-')
% %     
% %     hold on
% %     loglog(tau,y)
% %     title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
% %     xlabel('\tau (ms)')
% %     ylabel('Mean Square Distance (electrode)')
% %     
% %     str = {'p = ',num2str(p(1))};
% %     text(max(tau)+5,max(y)+5,str)
% % 
% %     pause
% %     close all
% end
% figure;
%     [MDsort,~,MDIdx] = unique(fullMD(:,2)) ;
%     
%     for idx = 1:length(MDsort)
%         meanMD(idx) = mean(fullMD(find(MDIdx==idx),1)) ;
%         stdMD(idx) = std(fullMD(find(MDIdx==idx),1)) ;
%     end
%     
%     % loglog(MSDsort,meanMSD,'.-')
%     errorbar(MDsort,meanMD,stdMD)
%     set(gca,'yscale','log','xscale','log')
%     hold on
%     
%     p = polyfit(log(MDsort(1:20)'),log(meanMD(1:20)),1) ;
%     y = exp(polyval(p,log(tau))) ;
%     slope = p(1) ;
%     loglog(tau,y)
%     title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
%     xlabel('\tau (ms)')
%     ylabel('Mean Square Distance (electrode)')
%     
%     str = {'p = ',num2str(p(1))};
%     text(max(tau)+5,max(y)+5,str)
% 
% figure
% hist(slopeAll,20)
% title('distribution of MSD for my144')
% xlabel('MSD')

%% displacement distribution
close all
center = WCentroids ;
deltaTime =  [2,4,8,16,32]; %,100] ;
iTa = 0.45 ;   % 0.26
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    DisplaceNorm = [] ;
for iBurst = 1:size(center,2)
    posCenter = center{iBurst} ;
    DisplaceTemp = sqrt(sum((posCenter(deltaT+1 :end,:) - posCenter(1: end-deltaT,:)).^2,2)) ;
    DisplaceTempNorm = DisplaceTemp/(deltaT^iTa) ;
    Displace = [Displace;DisplaceTemp] ;
    DisplaceNorm = [DisplaceNorm;DisplaceTempNorm] ;
end
DisplaceNormTemp = DisplaceNorm ;
%DisplaceNormTemp(DisplaceNormTemp<0.2) = [] ;
%DisplaceNormTemp(DisplaceNormTemp>20) = [] ;
[n,x] = hist(DisplaceNormTemp,200) ;
loglog(x,n/sum(n),'o')
hold on
legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
end
 legend(legendInfo)
 
% fit Gaussian
% pd = fitdist(Displace,'stable')
pd = fitdist(DisplaceNormTemp,'Normal')
y = pdf(pd,x) ;
hold on ; loglog(x,y/100,'r--');
legend([legendInfo,'Gaussian'])
xlabel('displacement (electrode)')
ylabel('probability')
title('Displacement of different time increment')

%% velocity autocorrelation

