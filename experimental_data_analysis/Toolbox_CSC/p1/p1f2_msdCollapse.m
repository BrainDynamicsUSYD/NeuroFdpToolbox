function p1f2_msdCollapse(WCentroids,superIdx)
%%
% figure;
center = WCentroids ;
deltaTime = [2,4,8,16,32] ; %,100]; %30 ;   
% iTa = 0.40 ;
iTa = 0.60 ;

legendInfo = [] ;
count = 0 ;
for delta = 1:length(deltaTime) %:-1:1
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    DisplaceNorm = [] ;
for iBurst = 1:length(superIdx) %  size(center,2)  % 
    curBurst = superIdx(iBurst) ;
    % curBurst = iBurst ;
    posCenter = center{curBurst} ;
    % DisplaceTemp = sqrt(sum((posCenter(deltaT+1 :end,:) - posCenter(1: end-deltaT,:)).^2,2)) ;
    DisplaceTemp = (posCenter(deltaT+1 :end,1) - posCenter(1: end-deltaT,1) ) ;
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
DisplaceNorm(DisplaceNorm==0) = [] ;
numPts = 200 ;
% nEdge = logspace(log10(min(DisplaceNorm)),log10(max(DisplaceNorm)),numPts) ;

[n,x] = histcounts(DisplaceNorm,numPts,'normalization','pdf') ;
% loglog(x(2:end),n,'o')
semilogy(0.5*(x(1:end-1)+x(2:end)),n,'o')
hold on
legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
end

fitVar = allDisplaceTemp ; % allDisplaceTemp ;  
pd = fitdist(fitVar,'stable') 
x = linspace(x(1)-1,x(end)+2,length(n)*2) ;
y = pdf(pd,x) ;
hold on ; % loglog(x,y,'r-');
semilogy(x,y,'r-','lineWidth',2);

pd2 = fitdist(fitVar,'normal') 
y = pdf(pd2,x) ;
hold on ; % loglog(x,y,'k--');
semilogy(x,y,'k--','lineWidth',1);

legend('t = 30 ms','Gaussian')
legendInfo{delta+1} = ['\alpha stable fit'] ;
legendInfo{delta+2} = ['Gaussian fit'] ;
legend(legendInfo)
xlabel('\rho (scaled displacement)')
ylabel('probability')
% title('Displacement distribution')

