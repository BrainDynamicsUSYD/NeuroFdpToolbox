function stableFit(sigDist,numPts)

sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
%paraDist = fitdist(sigDist','Stable') 
nEdge = logspace(log10(min(sigDist)),log10(max(sigDist)),numPts) ;

[x,n] = histcounts(sigDist,nEdge,'normalization','pdf') ;
n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ;
y0 = 2*pdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
% figure;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
loglog(n0,x,'b.','MarkerSize',12)
hold on
loglog(n0,y0,'r','lineWidth',2)
legend('original data','\alpha stable fit','Location','southwest')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

[R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable(sigDist) ;