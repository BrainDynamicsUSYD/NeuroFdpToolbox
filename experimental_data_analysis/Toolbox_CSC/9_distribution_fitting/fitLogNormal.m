close all
[PARMHAT,b] = lognfit(Duration) 

[n,x] = histcounts(Duration,400,'normalization','pdf') ;
semilogx(x(2:end),n,'.','MarkerSize',12)
hold on
xNew = linspace(x(2),x(end),1000) ;
y = lognpdf(xNew,PARMHAT(1),PARMHAT(2)) ;
semilogx(xNew,y,'r')
title('Duration distribution')
xlabel('ms')
ylabel('probability')
legend('orignal data','log-normal fit')

%% alpha stable
figure
pd = fitdist(Duration','stable') ;

[n,x] = histcounts(Duration,800,'normalization','pdf') ;
% [n,x] = hist(Duration,400,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
% semilogy(x,n/(sum(n)),'.','MarkerSize',12)
hold on
% xNew = linspace((x(2)+x(1))/2,(x(end-1)+x(end))/2,1000) ;
xNew = linspace((x(1)),(x(end)),1000) ;
y = pdf(pd,xNew) ;
semilogy(xNew,y,'r')
title('Duration distribution')
xlabel('ms')
ylabel('probability')

legend('original data','\alpha stable distribution fit')

%% exponential stable
figure
pd = fitdist(Duration','exp') ;

% [n,x] = histcounts(Duration,800,'normalization','pdf') ;
[n,x] = hist(Duration,100,'normalization','pdf') ;
% loglog((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
plot(x,n/(sum(n)),'.','MarkerSize',12)
hold on
% xNew = linspace((x(2)+x(1))/2,(x(end-1)+x(end))/2,1000) ;
xNew = linspace((x(1)),(x(end)),1000) ;
y = pdf(pd,xNew) ;
plot(xNew,y,'r')
title('Duration distribution')
xlabel('ms')
ylabel('probability')

legend('original data','exponential distribution fit')

%% log-normal
figure
[PARMHAT,b] = lognfit(patternScale) ;

[n,x] = histcounts(patternScale,400,'normalization','pdf') ;
semilogx(x(2:end),n,'.','MarkerSize',12)
hold on
xNew = linspace(x(2),x(end),1000) ;
y = lognpdf(xNew,PARMHAT(1),PARMHAT(2)) ;
semilogx(xNew,y,'r')
title('Scale distribution')
xlabel('electrode')
ylabel('probability')

legend('original data','log-normal fit')

%% alpha stable
figure
pd = fitdist(patternScale','stable') ;

[n,x] = histcounts(patternScale,400,'normalization','pdf') ;
loglog(x(2:end),n,'.','MarkerSize',12)
hold on
xNew = linspace(x(2),x(end),1000) ;
y = pdf(pd,xNew) ;
loglog(xNew,y,'r')
title('Scale distribution')
xlabel('electrode')
ylabel('probability')

legend('original data','\alpha stable distribution fit')

%%
figure
centIntervalTemp = centInterval ;
centIntervalTemp(centInterval<=0) = [];
PARMHAT = lognfit(centIntervalTemp) ;

[n,x] = histcounts(centIntervalTemp,400,'normalization','pdf') ;
semilogx((x(2:end)+x(1:end-1))/2,n,'.','markersize',12)
hold on
y = lognpdf(x,PARMHAT(1),PARMHAT(2)) ;
semilogx(x,y)
title('Interval distribution')
xlabel('time(ms)')
ylabel('probability')

legend('original data','log-normal fit')

%% alpha stable
close all
figure
centIntervalTemp = centInterval ;
centIntervalTemp(centInterval<0) = [] ;
sigPro = centIntervalTemp ;
pd = fitdist(sigPro','stable') 

[n,x] = histcounts(sigPro,100,'normalization','pdf') ;
semilogy((x(1:end-1)+x(2:end))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace(x(2),x(end),1000) ;
y = pdf(pd,xNew) ;
semilogy(xNew,y,'r')
title('Interval distribution')
xlabel('ms')
ylabel('probability')

% hold on
% pd = fitdist(sigPro','halfnormal') 
% y = pdf(pd,xNew) ;
% semilogy(xNew,y,'k--')

legend('original data','\alpha stable distribution fit')

%% exponential
close all
figure
centIntervalTemp = centInterval ;
centIntervalTemp(centInterval<0) = [] ;
sigPro = centIntervalTemp ;
pd = fitdist(sigPro','exp') 

[n,x] = histcounts(sigPro,100,'normalization','pdf') ;
plot((x(1:end-1)+x(2:end))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace(x(2),x(end),1000) ;
y = pdf(pd,xNew) ;
plot(xNew,y,'r')
title('Interval distribution')
xlabel('ms')
ylabel('probability')

% hold on
% pd = fitdist(sigPro','halfnormal') 
% y = pdf(pd,xNew) ;
% semilogy(xNew,y,'k--')

legend('original data','exponential')
