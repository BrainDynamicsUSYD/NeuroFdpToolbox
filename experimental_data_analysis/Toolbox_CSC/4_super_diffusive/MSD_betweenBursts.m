% this function is for characterizing the flights between bursts
%
% author: Xian Long    supervisor: Pulin Gong
%
%%
load('FullBurst95%ma027_032_StartminTime30msJune08_11:21')

%% histogram of flight time
[n,x] = hist(abs(centInterval),400) ;
loglog(x,n)

%% fluctuation of displacement
flucY = [] ;
tau = [] ;
countT = 1 ;
for iTau = 1:max(centInterval)
    idx = find(centInterval == iTau) ;
    if length(idx)<10
        continue
    else
        deltaY = distCent(idx) ;
        flucY(countT) = sqrt(mean(deltaY.^2)-mean(deltaY).^2) ;
        tau(countT) = iTau ;
        countT = countT+1 ;
        
    end
    
end
%%
% plot(tau,flucY)
loglog(tau,flucY)

%% check the burst coverage
burstArea = zeros(1,max(rangeFrame(:,2))) ;
for iBurst = 1:size(rangeFrame,1)
    tempRange = (rangeFrame(iBurst,1):rangeFrame(iBurst,2)) ;
    burstArea(tempRange) = burstArea(tempRange) +1 ;
    
    
    
end
plot(burstArea)
length(find(burstArea == 0))
length(find(burstArea == 1))
length(find(burstArea == 2))
length(find(burstArea == 3))
length(find(burstArea == 4))
length(find(burstArea == 5))
length(find(burstArea == 6))


%%


