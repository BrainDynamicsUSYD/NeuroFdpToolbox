%% This functions visualise the basic statistics after detecting the burst 
% patterns.
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis[n, xout] = hist(Duration,200);
[n, xout] = hist(Duration,200);
[r,p] = corrcoef([Duration(1:end-1)',peakAmp(1:end-1)',sumAmp(1:end-1)',...
    patternScale(1:end-1)',centInterval(2:end)',distCent(2:end)']) 
%%
close all
powerIndx = 1:0.05:5 ;
x = 10.^powerIndx ;

[n, xout] = hist(patternScale,x);
%[n, xout] = hist(patternScale,1000);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Size distribution over 300s')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(25:end)),log(nFit(25:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
    
figure
powerIndx = 1:0.05:3 ;
x = 10.^powerIndx ;
[n, xout] = hist(Duration,x);
% [n, xout] = hist(Duration,1000) ;

n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Duration distribution over 300s')
xlabel('Lifetime of a burst (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout),log(nFit),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
    

%%
close all
clear
load('FullBurstma025_03_StartminTime30msMay08_22:25.mat')
CentroidsNew{1} = Centroids ;
WCentroidsNew{1} = WCentroids ;
distCrossBurst{1} = distCent ;

load('FullBurstma027_032_StartminTime30msMay08_19:18.mat')
CentroidsNew{2} = Centroids ;
WCentroidsNew{2} = WCentroids ;
distCrossBurst{2} = distCent ;

load('FullBurstmy144_101_StartminTime30msMay08_18:02.mat')
CentroidsNew{3} = Centroids ;
WCentroidsNew{3} = WCentroids ;
distCrossBurst{3} = distCent ;

load('FullBurstmy147_53_StartminTime30msMay08_21:43.mat')
CentroidsNew{4} = Centroids ;
WCentroidsNew{4} = WCentroids ;
distCrossBurst{4} = distCent ;

%%
load('FullBurstSur1ma025_03_StartminTime30msMay09_07:34.mat')
SCentroidsNew{1} = Centroids ;
SWCentroidsNew{1} = WCentroids ;
SdistCrossBurst{1} = distCent ;

load('FullBurstSur1ma027_032_StartminTime30msMay09_05:44.mat')
SCentroidsNew{2} = Centroids ;
SWCentroidsNew{2} = WCentroids ;
SdistCrossBurst{2} = distCent ;

load('FullBurstSur1my144_101_StartminTime30msMay09_08:12.mat')
SCentroidsNew{3} = Centroids ;
SWCentroidsNew{3} = WCentroids ;
SdistCrossBurst{3} = distCent ;

load('FullBurstSur1my147_53_StartminTime30msMay09_05:57.mat')
SCentroidsNew{4} = Centroids ;
SWCentroidsNew{4} = WCentroids ;
SdistCrossBurst{4} = distCent ;

%%
figure;
for i = 1:4
    
center = WCentroidsNew{i} ;
distCentBurstSum = [] ;
xdistCentBurstSum = [] ;
ydistCentBurstSum = [] ;
for iBurst = 1:size(center,2)
    centreBurst = center(center(:,iBurst,1)~=0,iBurst,:) ;
    centreBurst = squeeze(centreBurst) ;
    % distCentBurst = sqrt( sum( ( centreBurst(1:end-1,:)-centreBurst(2:end,:) ).^2,2) ) ;
    distCentBurst = sqrt( sum( (diff(centreBurst)).^2,2 ) ) ;
    distCentBurstSum = [ distCentBurstSum;distCentBurst] ;
    
    xdistCentBurst =  centreBurst(1:end-1,1)-centreBurst(2:end,1) ;
    xdistCentBurstSum = [ xdistCentBurstSum;xdistCentBurst] ;
    
    ydistCentBurst =  centreBurst(1:end-1,2)-centreBurst(2:end,2) ;
    ydistCentBurstSum = [ ydistCentBurstSum;ydistCentBurst] ;
end
% figure
% subplot(3,1,1)
% hist(distCentBurstSum,4000)
% title('Distribution of distance of centre')
% xlim([-2,2])
% ylim([0,4000])
% subplot(3,1,2)
% hist(xdistCentBurstSum,2000)
% title('Distribution of variation at x axis')
% xlim([-2,2])
% ylim([0,4000])
% subplot(3,1,3)
% hist(ydistCentBurstSum,2000)
% title('Distribution of variation at y axis')
% xlim([-2,2])
% ylim([0,4000])
% xlabel('electrode')

% comment this line to plot distribution of only within burst
distCentBurstSum = [distCentBurstSum;distCrossBurst{i}' ];

subplot(4,2,2*i-1)
[n,xout] = hist(distCentBurstSum,100) ;
plot(xout,n,'.','MarkerSize',12)
% title('Distribution of distance of centre')
%xlim([0 14])
ylim([0 1e4])

subplot(4,2,2*i)
[n,xout] = hist(distCentBurstSum,100) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
% title('Distribution of distance of centre')
% xlim([0 14])
ylim([0 1e4])

% % fit a levy distribution
% pd = fitdist(distCentBurstSum,'stable') ;
% x_values = -4:0.001:4;
% y = pdf(pd,x_values);
% plot(x_values,y,'r','LineWidth',2)

end

%%

DurationAll = Duration ;
patternScaleAll = patternScale ;
surScale = zeros(size(patternScaleAll,2),1) ;

%%
for nSur = 1:size(patternScaleAll,2)
    surScale(nSur) = mean(patternScaleAll(patternScaleAll(:,nSur)~=0,nSur) ) ;
       
end
hist(surScale,20)
meanScale = mean(patternScale) ;


%% hist of interval
figure;

hist(centInterval,80)
title('distribution of time interval between bursts')
xlabel('time (ms)')
xlim([min(centInterval) max(centInterval)])

%%
figure
powerIndx = 1:0.05:5 ;
x = 10.^powerIndx ;

[n,xout] = hist((centInterval(centInterval>0)),x) ;
% [n,xout] = hist(abs(centInterval(centInterval>0)),50) ;
% [n,xout] = hist(abs(centInterval),100) ;
% temp = centInterval ;
% temp(temp<0) = 0 ;
% [n,xout] = hist(temp,100) ;

loglog(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')

nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
nFit(isnan(log(n))) = [] ;
xout(isnan(log(n))) = [] ;
p = polyfit(log(xout(15:end)),log(nFit(15:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y)-1,str) 

title('distribution of time interval between bursts')
xlabel('time (ms)')


%% hist of interval
figure;
hist(distCent,10)
title('distribution of jumping distance between bursts')
xlabel('time (ms)')
xlim([min(distCent) max(distCent)])

% figure
% [n,xout] = hist(abs(centInterval(centInterval>0)),100) ;
% % [n,xout] = hist(abs(centInterval),100) ;
% % temp = centInterval ;
% % temp(temp<0) = 0 ;
% % [n,xout] = hist(temp,100) ;
% 
% loglog(xout,n,'o')
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
% nFit(isnan(log(n))) = [] ;
% xout(isnan(log(n))) = [] ;
% p = polyfit(log(xout(4:35)),log(nFit(4:35)),1) ;
% y = exp(polyval(p,log(xout))) ;
% hold on
% loglog(xout,y)
% str = {'p = ',num2str(p(1))};
% text(max(xout)-1,max(y)-1,str) 
% 
% title('distribution of time interval between bursts')
% xlabel('time (ms)')

%% where the burst is in time
close all
middlePoint = rangeFrame(:,1)' + Duration/2 ;
rangePoint = Duration ;
x = 1:length(middlePoint) ;
errorbar(middlePoint,x,rangePoint,'horizontal','o')
xlabel('time (ms)')
ylabel('burst')
title('Temporal location of bursts')
xlim([0e3+1,10e3])


%% histogram of displacement between bursts
close all
numHistPts = 40 ;
hist(distCent,numHistPts)
title('Histogram of displacement between bursts')
xlabel('electrode')
ylabel('count')
[n,x] = hist(distCent,numHistPts) ;
figure;
loglog(x,n,'o')
figure;
semilogx(x,n)

%% histogram of displacement between burst centre
center = WCentroids ;
meanLoc = zeros(size(center,2),2) ;
for iBurst = 1:size(center,2)
    meanLoc(iBurst,:) = mean(WCentroids{iBurst},1) ;
end

distMeanLoc = sqrt(sum((diff(meanLoc,[],1).^2),2)) ;
figure;
hist(distMeanLoc,numHistPts)
[n2,x2] = hist(distMeanLoc,numHistPts) ;
figure;
loglog(x2,n2,'o')
figure;
semilogx(x2,n2)

%% histogram of burst location
figure;
hist(meanLoc(:,1),20)
figure;
hist(meanLoc(:,2),20)
