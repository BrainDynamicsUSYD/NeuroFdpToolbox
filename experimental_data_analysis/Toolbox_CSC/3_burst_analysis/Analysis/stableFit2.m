function stableFit2(sigDist,numPts)
sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
%paraDist = fitdist(sigDist','Stable') 
[x,n] = histcounts(sigDist,50,'normalization','pdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = 2*pdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
semilogy(n0,x,'b.','MarkerSize',12)
hold on
semilogy(n0,y0,'r','lineWidth',2)
legend('original data','\alpha stable fit')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

[R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable(sigDist) ;