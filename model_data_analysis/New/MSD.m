
%% check how many bursts at a same time

burst = zeros(size(Centroids,1),1) ;

for iBurst = 1:size(Centroids,2)
    temp = zeros(size(Centroids,1),1) ;
    temp(Centroids(:,iBurst,1)~=0)  = 1 ;
    burst = burst + temp ;
    
end
plot(burst)
%% select only 1 burst
center = WCentroids ;
% slope = [] ;
for iBurst = 20
    Trajectory = [] ;
    temp = zeros(size(center,1),1) ;
    TrajectoryTemp = squeeze(center(center(:,iBurst,1)~=0,iBurst,:)) ;
    TrajectoryTemp(:,3) = find(center(:,iBurst,1)~=0) ;
    Trajectory = [Trajectory;TrajectoryTemp] ;
    

    [MSD,tau] = get_MSD(Trajectory(1:end,:)) ;%,grid_size);
     loglog(tau,MSD,'.-')
    p = polyfit(log(tau(1:10)),log(MSD(1:10)),1) ;
    y = exp(polyval(p,log(tau))) ;
     hold on
    loglog(tau,y)
    title('Mean Square Distance versus \tau within a burst')
    xlabel('\tau (ms)')
    ylabel('Mean Square Distance (electrode)')
    
    str = {'p = ',num2str(p(1))};
    text(max(tau)+5,max(y)+5,str)
    
    slope(iBurst) = p(1) ;
end
hist(slope)
