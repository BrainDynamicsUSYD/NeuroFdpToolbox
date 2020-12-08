%% levy walk simulation
mu = 1.15 ;
X = rand(10000,1) * pi - pi/2 ;
X_prime = rand(10000,1) ;
Y = -log(X_prime) ;
Z = sin( (mu-1)*X )./(cos(X).^(mu-1)) .* (cos((2-mu)*X)./Y).^((2-mu)/(mu-1)) ;

X0 = 5;
Y0 = 5;
traj = [] ;
for iTime = 1:10000 
    direction = rand*2*pi ;
    XNew = X0+Z(iTime)*cos(direction) ;
    YNew = Y0+Z(iTime)*sin(direction) ;
    X0 = XNew ;
    Y0 = YNew ;
    coord = [XNew,YNew,iTime] ;
    traj = [traj; coord];
    % plot(XNew,YNew,'o')
    % hold on
end

%% select only 1 burst
Trajectory = traj ;
% % select bursts of duration more than 100ms
% center = [] ;
% count = 1 ;
% for iBurst = 1:length(Duration)
%     if Duration(iBurst)>100
%         center{count} = WCentroids{iBurst} ;
%         count = count+1 ;
%     end
% end

slope = [] ;
fullMSD = [] ;
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 0 ;

for iBurst = 1% :size(Trajectory,1)    
    [MSD,tau] = get_MSD(Trajectory) ;%,grid_size);
    
    tempMSD = [MSD,tau] ;
    fullMSD = [fullMSD;tempMSD] ;

    
    pAll = polyfit(log(tau(1:20)),log(MSD(1:20)),1) ;
    y = exp(polyval(pAll,log(tau))) ;
    slopeAll(iBurst) = pAll(1) ;
    
    if pAll(1)>1.1
        msdSuper = [msdSuper;tempMSD] ; 
        countSuper = countSuper + 1;
    end
    %plot
%     loglog(tau,MSD,'.-')
%     
%     hold on
%     loglog(tau,y)
%     title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
%     xlabel('\tau (ms)')
%     ylabel('Mean Square Distance (electrode)')
%     
%     str = {'p = ',num2str(p(1))};
%     text(max(tau)+5,max(y)+5,str)
% 
%     pause
%     close all
end
figure;
    [MSDsort,~,MSDIdx] = unique(fullMSD(:,2)) ;
    
    for idx = 1:length(MSDsort)
        meanMSD(idx) = mean(fullMSD(find(MSDIdx==idx),1)) ;
        stdMSD(idx) = std(fullMSD(find(MSDIdx==idx),1)) ;
    end
    
    % loglog(MSDsort,meanMSD,'.-')
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
hist(slopeAll,20)
title('distribution of MSD for my144')
xlabel('MSD')


% [n,x] = hist(X,1000) ;
%loglog(x,n,'o') ;
