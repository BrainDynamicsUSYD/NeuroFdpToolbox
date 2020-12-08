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


%% velocity autocorrelation

center = WCentroids ;
deltaTime = 0:120 ; % [2,6,12,25,50,100] ;
Vall = zeros(1,length(deltaTime)) ;
dT = 1 ;
count = 0 ;

for iBurst = 1:size(center,2)
    posCenter = center{iBurst} ;
    Velocity = (posCenter(dT+1 :end,:) - posCenter(1: end-dT,:))/dT ;
    for delta = 1:length(deltaTime)
        deltaT = deltaTime(delta) ;
        V_deltaT(delta) = mean(sum(Velocity(deltaT+1:end,:).*Velocity(1:end-deltaT,:),2 )) ;
        if V_deltaT(1) == 0
            V_deltaT(2:end) = zeros(1,length(deltaTime)-1) ;
            V_deltaT(1) = 1 ;
            count = count-1 ;
            break
        end
    end
    V_normalised = V_deltaT(1:end)/V_deltaT(1) ;

    loglog(V_normalised,'o')
title(['burst ',num2str(iBurst)])
     pause
     close all
% Vall = Vall+ V_normalised ;
% count = count+1 ;

end
% loglog(Vall/count,'o')
%% displacement autocorrelation
center = WCentroids ;
acfAll = [] ;
KOld = zeros(20,1) ;
for iBurst = 1:size(center,2)
    posCenter = center{iBurst} ;
    K= [] ;
    if size(posCenter,1)<36
        continue
    else
    for tau = 1:30
        deltaT = 1:5 ;
            DisplaceTemp1 = posCenter(deltaT+1 ,:) - posCenter(1,:) ;
            DisplaceTemp2 = posCenter(deltaT+1+tau ,:) - posCenter(1+tau,:) ;
            K(tau) = mean(sum(DisplaceTemp1.*DisplaceTemp2,2)) ;
        
    end
    end
   loglog(K,'o')    
KFull = K+KOld ; 

KOld = K ; 

%    autocorr(DisplaceTemp)
%     xVelocity = (posCenter(deltaT+1 :end,1) - posCenter(1: end-deltaT,1))/deltaTime ;
%     yVelocity = (posCenter(deltaT+1 :end,2) - posCenter(1: end-deltaT,2))/deltaTime ;
%     autocorr(xVelocity) 
%      autocorr(xVelocity) 
    %acf = abs(acf1+acf2) ;
%     if length(acf)>length(acfAll)
%         acfAll = [acfAll; zeros(length(acf)-length(acfAll),1) ] ;
%     elseif length(acf)<length(acfAll)
%         acf = [acf; zeros(length(acfAll)-length(acf) ,1) ]
%     end
    %acfAll = acfAll+acf ;
%     loglog(acf)
%     figure;
%     loglog(acf1)
%     figure;
%     loglog(acf2)
    pause
%    close all
    % Displace = [Displace;DisplaceTemp] ;
end