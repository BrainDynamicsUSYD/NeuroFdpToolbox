Trajectory = [] ;
for iBurst = 1:size(rangeFrame,1)  
    tempTrajectory = [] ;
    center = WCentroids{iBurst} ;

    tempTrajectory(:,1:2) = center ;
    tempTrajectory(:,3) = rangeFrame(iBurst,1):rangeFrame(iBurst,2) ;
    Trajectory = [Trajectory ; tempTrajectory] ;
end
% [MSD,tau] = get_MSD(Trajectory(0*fsTemporal+1:fix(10*fsTemporal),:)) ;%,grid_size);
[MSD,tau] = get_MSD(Trajectory) ;
% loglog(tau,MSD)
close all
%%
loglog(tau,MSD,'o-')

idxStart = find(tau==2) ;
idxEnd = find(tau==22) ;    %74

        fitRange = (idxStart:idxEnd) ;
        [pAll,S] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
        y = exp(polyval(pAll,log(tau(fitRange)))) ;
        p_temp = pAll(1)
        hold on 
        loglog(tau(fitRange),y,'r')
        title(['MSD for 30-40s p=',num2str(p_temp)]) 
        xlabel('tau')
        ylabel('electrode')

%% only between bursts
Trajectory = [] ;
for iBurst = 1:size(rangeFrame,1)  
    tempTrajectory = [] ;
    center = WCentroids{iBurst} ;

    tempTrajectory(:,1:2) = center(end,:) ;
    tempTrajectory(:,3) = rangeFrame(iBurst,2) ;
    Trajectory = [Trajectory ; tempTrajectory] ;
end

% [MSD,tau] = get_MSD(Trajectory(0*fsTemporal+1:fix(10*fsTemporal),:)) ;%,grid_size);
[MSD,tau] = get_MSD(Trajectory) ;
% loglog(tau,MSD)
close all
%
loglog(tau,MSD,'o-')

idxStart = find(tau==2) ;
idxEnd = find(tau==22) ;    %74

fitRange = (idxStart:idxEnd) ;
[pAll,S] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
y = exp(polyval(pAll,log(tau(fitRange)))) ;
p_temp = pAll(1)
hold on
loglog(tau(fitRange),y,'r')
title(['MSD for 30-40s p=',num2str(p_temp)])
xlabel('tau')
ylabel('electrode')
