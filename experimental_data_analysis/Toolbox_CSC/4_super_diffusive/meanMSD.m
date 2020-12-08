center = WCentroids ;

for iBurst = 1:size(center,2)
    tempLoc = center{iBurst} ;
    distLoc(iBurst) = mean(sqrt(sum(diff(tempLoc).^2,2))) ;
end
%%
hist(distLoc,20)

figure;
distLoc(distLoc==0) = [] ;
PARMHAT = lognfit(distLoc) ;

[n,x] = hist(distLoc,40)
semilogx(x,n/sum(n),'o')
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y/sum(y))
title('Duration')
xlabel('ms')
ylabel('probability')