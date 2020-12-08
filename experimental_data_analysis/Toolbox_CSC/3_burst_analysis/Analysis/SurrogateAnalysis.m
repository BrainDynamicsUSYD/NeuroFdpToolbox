%% This function analyzes the results (statistics) of surrogate data
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. raw data
% load('FullBurst95%ma027_032_StartminTime30msMay24_21:18.mat')
numBurst_raw = size(Duration,2) ;
duration_raw = Duration ;
scale_raw = patternScale ;
interval_raw = centInterval ;
distCent_raw = distCent ;

%% 2. amplitude according to Gaussian (AAFT)
% load('FullBurst95%SurAAFTma027_032_StartminTime30msMay28_18:58.mat')
numBurst_AAFT = size(Duration,2) ;
duration_AAFT = Duration ;
scale_AAFT = patternScale ;
interval_AAFT = centInterval ;
distCent_AAFT = distCent ;

%% 3. spatial randomized
load('SpaceOnlyFullBurst95%Surma027_032_StartminTime30msMay29_09:03.mat')
numBurst_space = size(Duration,2) ;
duration_space = Duration ;
scale_space = patternScale ;
interval_space = centInterval ;
distCent_space = distCent ;

%% 4. phase randomized
% load('FullBurst95%Surma027_032_StartminTime30msMay28_17:18.mat')
numBurst_phase = size(Duration,2) ;
duration_phase = Duration ;
scale_phase = patternScale ;
interval_phase = centInterval ;
distCent_phase = distCent ;

%% 5. phase and spatial randomized
% load('SpaceFullBurst95%Surma027_032_StartminTime30msMay28_18:26.mat')
numBurst_phase_space = size(Duration,2) ;
duration_phase_space = Duration ;
scale_phase_space = patternScale ;
interval_phase_space = centInterval ;
distCent_phase_space = distCent ;

%% 6. Gaussian
% load('SpaceOnlyFullBurst95%SurGaussianNoise_StartminTime30msMay29_09:17.mat')
numBurst_Gaussian = size(Duration,2) ;
duration_Gaussian = Duration ;
scale_Gaussian = patternScale ;
interval_Gaussian = centInterval ;
distCent_Gaussian = distCent ;


%% scale
subplot(3,2,1)
[n, xout] = hist(scale_raw,numBurst_raw);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Experimental data')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])    


subplot(3,2,2)
[n, xout] = hist(scale_AAFT,numBurst_AAFT);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Amplitude adjusted Fourier transform')
xlabel('Total sites in a burst')
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
xlim([1 1e5])    


subplot(3,2,3)
[n, xout] = hist(scale_phase,numBurst_phase);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized phase')
xlabel('Total sites in a burst')
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
xlim([1 1e5])    
   


subplot(3,2,4)
[n, xout] = hist(scale_space,numBurst_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])    


subplot(3,2,5)
[n, xout] = hist(scale_phase_space,numBurst_phase_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space and phase')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])  

subplot(3,2,6)
[n, xout] = hist(scale_Gaussian,numBurst_Gaussian);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Gaussian noise')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5]) 

%% interval

subplot(3,2,1)
interval_raw_temp = interval_raw ;
interval_raw_temp(interval_raw<0) = 0 ;
[n, xout] = hist(interval_raw_temp,numBurst_raw);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Experimental data')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
p = polyfit( log(xout),log(nFit) ,1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])    
ylim([1e-4 1e-1])


subplot(3,2,2)
interval_AAFT_temp = interval_AAFT ;
interval_AAFT_temp(interval_AAFT<0) = 0 ;
[n, xout] = hist(interval_AAFT_temp,numBurst_AAFT);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Amplitude adjusted Fourier transform')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout),log(nFit),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])    
ylim([1e-4 1e-1])


subplot(3,2,3)
interval_phase_temp = interval_phase ;
interval_phase_temp(interval_phase<0) = 0 ;
[n, xout] = hist(interval_phase_temp,numBurst_phase);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized phase')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout),log(nFit),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])    
ylim([1e-4 1e-1])
   


subplot(3,2,4)
interval_space_temp = interval_space ;
interval_space_temp(interval_space<0) = 0 ;
[n, xout] = hist(interval_space_temp,numBurst_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])    
ylim([1e-4 1e-1])


subplot(3,2,5)
interval_phase_space_temp = interval_phase_space ;
interval_phase_space_temp(interval_phase_space<0) = 0 ;
[n, xout] = hist(interval_phase_space_temp,numBurst_phase_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space and phase')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])  
ylim([1e-4 1e-1])


subplot(3,2,6)
interval_Gaussian_temp = interval_Gaussian ;
interval_Gaussian_temp(interval_Gaussian<0) = 0 ;
[n, xout] = hist(interval_Gaussian_temp,numBurst_Gaussian);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Gaussian noise')
xlabel('Interval (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout),1e-3,str)    
xlim([0 1e5])  
ylim([1e-4 1e-1])
  
%% duration

subplot(3,2,1)
[n, xout] = hist(duration_raw(1:400),400);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Experimental data')
xlabel('Duration (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])    


subplot(3,2,2)
[n, xout] = hist(duration_AAFT,numBurst_AAFT);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Amplitude adjusted Fourier transform')
xlabel('Duration (ms)')
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
xlim([1 1e5])    


subplot(3,2,3)
[n, xout] = hist(duration_phase,numBurst_phase);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized phase')
xlabel('Duration (ms)')
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
xlim([1 1e5])    
   


subplot(3,2,4)
[n, xout] = hist(duration_space,numBurst_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space')
xlabel('Duration (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])    


subplot(3,2,5)
[n, xout] = hist(duration_phase_space,numBurst_phase_space);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Randomized space and phase')
xlabel('Duration (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])  

subplot(3,2,6)
[n, xout] = hist(duration_Gaussian,numBurst_Gaussian);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Gaussian noise')
xlabel('Duration (ms)')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(1:end)),log(nFit(1:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([1 1e5])  


%%


[alpha, xmin, L]=plfit(scale_raw) ;

h=plplot(scale_raw, xmin, alpha)

%% scale
subplot(3,2,1)
[alpha, xmin, L]=plfit(scale_raw) ;

h=plplot(scale_raw, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])    
title('Experimental data')


subplot(3,2,2)
[alpha, xmin, L]=plfit(scale_AAFT) ;

h=plplot(scale_AAFT, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])  
title('Amplitude adjusted Fourier transform')


subplot(3,2,3)
[alpha, xmin, L]=plfit(scale_phase) ;

h=plplot(scale_phase, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5]) 
title('Randomized phase')



subplot(3,2,4)
[alpha, xmin, L]=plfit(scale_space) ;

h=plplot(scale_space, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space')


subplot(3,2,5)
[alpha, xmin, L]=plfit(scale_phase_space) ;

h=plplot(scale_phase_space, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space and phase')


subplot(3,2,6)
[n, xout] = hist(scale_Gaussian,numBurst_Gaussian);
[alpha, xmin, L]=plfit(scale_Gaussian) ;

h=plplot(scale_Gaussian, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Gaussian noise')

%% duration
subplot(3,2,1)
[alpha, xmin, L]=plfit(duration_raw) ;

h=plplot(duration_raw, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])    
title('Experimental data')


subplot(3,2,2)
[alpha, xmin, L]=plfit(duration_AAFT) ;

h=plplot(duration_AAFT, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])  
title('Amplitude adjusted Fourier transform')


subplot(3,2,3)
[alpha, xmin, L]=plfit(duration_phase) ;

h=plplot(duration_phase, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5]) 
title('Randomized phase')



subplot(3,2,4)
[alpha, xmin, L]=plfit(duration_space) ;

h=plplot(duration_space, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space')


subplot(3,2,5)
[alpha, xmin, L]=plfit(duration_phase_space) ;

h=plplot(duration_phase_space, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space and phase')


subplot(3,2,6)
[n, xout] = hist(duration_Gaussian,numBurst_Gaussian);
[alpha, xmin, L]=plfit(duration_Gaussian) ;

h=plplot(scale_Gaussian, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Gaussian noise')

%% interval
subplot(3,2,1)
interval_raw_temp = interval_raw(interval_raw>0)  ;
[alpha, xmin, L]=plfit(interval_raw_temp,'finite') ;

h=plplot(interval_raw_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])    
title('Experimental data')


subplot(3,2,2)
interval_AAFT_temp = interval_AAFT(interval_AAFT>0)  ;
[alpha, xmin, L]=plfit(interval_AAFT_temp,'finite') ;

h=plplot(interval_AAFT_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])  
title('Amplitude adjusted Fourier transform')


subplot(3,2,3)
interval_phase_temp = interval_phase(interval_phase>0) ;
[alpha, xmin, L]=plfit(interval_phase_temp,'finite') ;

h=plplot(interval_phase_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5]) 
title('Randomized phase')



subplot(3,2,4)
interval_space_temp = interval_space(interval_space>0) ;
[alpha, xmin, L]=plfit(interval_space_temp,'finite') ;

h=plplot(interval_space_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space')


subplot(3,2,5)
interval_phase_space_temp = interval_phase_space(interval_phase_space>0) ;
[alpha, xmin, L]=plfit(interval_phase_space_temp,'finite') ;

h=plplot(interval_phase_space_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Randomized space and phase')


subplot(3,2,6)
interval_Gaussian_temp = interval_Gaussian(interval_Gaussian>0) ;
[alpha, xmin, L]=plfit(interval_Gaussian_temp,'finite') ;

h=plplot(interval_Gaussian_temp, xmin, alpha)

str = {'p = ',num2str(alpha)};
text(xmin,1e-3,str)    
xlim([1 1e5])   
title('Gaussian noise')

%% Scale

[alpha, xmin, L]=plfit(scale_phase ,'finite')
h=plplot(scale_phase, xmin, alpha, 'ro')
str = {'slope = ',num2str(alpha)};
text(max(scale_phase),1e-3,str) 
legend('Surrogate data')
hold on

[alpha, xmin, L]=plfit(scale_raw,'finite')
h2=plplot(scale_raw, xmin, alpha,'bo')
% [p,gof]=plpva(patternScale, xmin)
str = {'slope = ',num2str(alpha)};
text(mean(scale_raw),1e-3,str)    
a=axes('position',get(gca,'position'),'visible','off');
legend(a,h2,'Raw data')
title('Size for ma027 (95% and 30ms)')

%% Duration
[alpha, xmin, L]=plfit(duration_phase ,'finite')
h=plplot(duration_phase, xmin, alpha, 'ro')
str = {'slope = ',num2str(alpha)};
text(max(duration_phase),1e-3,str) 
legend('Surrogate data')
hold on

[alpha, xmin, L]=plfit(Duration,'finite')
h2=plplot(Duration, xmin, alpha,'bo')
% [p,gof]=plpva(patternScale, xmin)
str = {'slope = ',num2str(alpha)};
text(mean(Duration),1e-3,str)    
a=axes('position',get(gca,'position'),'visible','off');
legend(a,h2,'Raw data')
title('Duration for ma027 (95% and 30ms)')

%% Interval
interval_phase(interval_phase<=0) = [] ;
[alpha, xmin, L]=plfit(interval_phase ,'finite')
h=plplot(interval_phase, xmin, alpha, 'ro') ;
str = {'slope = ',num2str(alpha)};
text(max(interval_phase),1e-3,str) 
legend('Surrogate data')
hold on

interval_phase(interval_phase<=0) = [] ;
[alpha, xmin, L]=plfit(interval_raw,'finite')
h2=plplot(interval_raw, xmin, alpha,'bo') ;
% [p,gof]=plpva(patternScale, xmin)
str = {'slope = ',num2str(alpha)};
text(max(interval_raw)-1000,1e-3,str)    
a=axes('position',get(gca,'position'),'visible','off');
legend(a,h2,'Raw data')
title('Interval for ma027 (95% and 30ms)')


%% coef

corrcoef(scale_raw(1:end-1),interval_raw(2:end))
corrcoef(scale_AAFT(1:end-1),interval_AAFT(2:end))
corrcoef(scale_phase(1:end-1),interval_phase(2:end))
corrcoef(scale_space(1:end-1),interval_space(2:end))

%%
% patternScaleOld = [] ;
% DurationOld = [] ;
% centIntervalOld = [] ;

patternScaleOld = [patternScaleOld,patternScale] ;

DurationOld = [DurationOld,Duration] ;

centIntervalOld = [centIntervalOld,patternScale] ;

%%
z = linspace(1e1,1e5,1000);
x = 0.5*(z(1:end-1)+z(2:end));
figure;
subplot(2,1,1)
[n, xout] = histcounts(scale_raw,z,'normalization','pdf');
% n=n/sum(n) ;
plot(x,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Experimental data')
xlabel('Total sites in a burst')
ylabel('Probability')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(10:2/3*length(nFit))),log(nFit(10:2/3*length(nFit))),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)    
xlim([30 5e5])   
ylim([1e-6 1e-2])
legend('original data','power law fitting')

subplot(2,1,2)
[n, xout] = histcounts(scale_phase,z,'normalization','pdf');
% n=n/sum(n) ;
plot(x,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Surrogate data')
xlabel('Total sites in a burst')
ylabel('Probability')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
% p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
% y = exp(polyval(p,log(xout))) ;
% hold on
% loglog(xout,y)
% str = {'p = ',num2str(p(1))};
% text(max(xout)-1,max(y),str)
% pd = fitdist(n','lognormal') ;
% nFit = pdf(pd,x) ;
% loglog(x,nFit) ;
hold on

% MLE of power law with cut off
xmin = 150 ;
sigIn = scale_phase ;
sigIn = sigIn(sigIn>xmin) ;
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;

c = histcounts(sigIn,linspace(1e1,1e5,1000),'Normalization','pdf');
% loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));


xlim([30 5e5])
ylim([1e-6 1e-2])
legend('original data','exponential fitting')


%%
z = linspace(20,1000,500);
x = 0.5*(z(1:end-1)+z(2:end));

figure;
subplot(2,1,1)
% [n, xout] = histcounts(duration_raw,z,'normalization','pdf');
% % n=n/sum(n) ;
% plot(x,n,'.','MarkerSize',12)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('Experimental data')
% xlabel('Duration (ms)')
% ylabel('Probability')
% nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
% p = polyfit(log(xout(10:2/3*length(nFit))),log(nFit(10:2/3*length(nFit))),1) ;
% y = exp(polyval(p,log(xout))) ;
% hold on
% loglog(xout,y)
% str = {'p = ',num2str(p(1))};
% text(max(xout)-1,max(y),str)    
% xlim([20 3e3])   
% ylim([1e-5 1])
% legend('original data','power law fitting')
% alpha stable fit

sigInProp = duration_raw ;
pd = fitdist(sigInProp','stable') ;

[n,x] = histcounts(sigInProp,400,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
semilogy(xNew,y,'r')
title('Duration distribution of experimental data')
xlabel('Time (ms)')
ylabel('probability')
xlim([20 7e3])
ylim([1e-6 1e-1])
legend('original data','\alpha stable distribution fit')

subplot(2,1,2)
% [n, xout] = histcounts(duration_phase,z,'normalization','pdf');
% % n=n/sum(n) ;
% plot(x,n,'.','MarkerSize',12)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('Surrogate data')
% xlabel('Duration (ms)')
% ylabel('Probability')
% nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
% p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
% y = exp(polyval(p,log(xout))) ;
% xlim([20 3e3])
% ylim([1e-5 1])
% % hold on
% % loglog(xout,y)
% % str = {'p = ',num2str(p(1))};
% % text(max(xout)-1,max(y),str)
% hold on
% 
% % MLE of power law with cut off
% xmin = 30 ;
% sigIn = duration_phase ;
% sigIn = sigIn(sigIn>xmin) ;
% pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
% [lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
% L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;
% 
% c = histcounts(sigIn,linspace(1e1,1e5,1000),'Normalization','pdf');
% % loglog(x,c,'.','markersize',20)
% hold on
% loglog(x,pf_exp(x, lambdaHat3));
% 
% xlim([20 3e3])
% ylim([1e-5 1])
% legend('original data','exponential fitting')

% alpha stable fit

sigInProp = duration_phase ;
% pd = fitdist(sigInProp','exponential') ;

[n,x] = histcounts(sigInProp,400,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
% semilogy(xNew,y,'r')
title('Duration distribution of randomized data')
xlabel('Time (ms)')
ylabel('probability')
xlim([20 6e3])   
ylim([1e-6 1e-1])

% legend('original data','\alpha stable distribution fit')


%%
z = linspace(20,1000,500);
x = 0.5*(z(1:end-1)+z(2:end));

interval_raw(interval_raw<=0) = [] ;
interval_raw(interval_raw<=0) = [] ;

figure;
subplot(2,1,1)
% [n, xout] = histcounts(duration_raw,z,'normalization','pdf');
% % n=n/sum(n) ;
% plot(x,n,'.','MarkerSize',12)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('Experimental data')
% xlabel('Duration (ms)')
% ylabel('Probability')
% nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
% p = polyfit(log(xout(10:2/3*length(nFit))),log(nFit(10:2/3*length(nFit))),1) ;
% y = exp(polyval(p,log(xout))) ;
% hold on
% loglog(xout,y)
% str = {'p = ',num2str(p(1))};
% text(max(xout)-1,max(y),str)    
% xlim([20 3e3])   
% ylim([1e-5 1])
% legend('original data','power law fitting')
% alpha stable fit

sigInProp = interval_raw ;
pd = fitdist(sigInProp','stable') ;

[n,x] = histcounts(sigInProp,400,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
semilogy(xNew,y,'r')
title('Interval distribution of experimental data')
xlabel('Time (ms)')
ylabel('probability')
xlim([20 4.5e3])   
ylim([1e-5 1e-1])
legend('original data','\alpha stable distribution fit')

subplot(2,1,2)
% [n, xout] = histcounts(duration_phase,z,'normalization','pdf');
% % n=n/sum(n) ;
% plot(x,n,'.','MarkerSize',12)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('Surrogate data')
% xlabel('Duration (ms)')
% ylabel('Probability')
% nFit = n ;
% nFit(~isfinite(log(n))) = [] ;
% xout(~isfinite(log(n))) = [] ;
% p = polyfit(log(xout(1:2/3*length(nFit))),log(nFit(1:2/3*length(nFit))),1) ;
% y = exp(polyval(p,log(xout))) ;
% xlim([20 3e3])
% ylim([1e-5 1])
% % hold on
% % loglog(xout,y)
% % str = {'p = ',num2str(p(1))};
% % text(max(xout)-1,max(y),str)
% hold on
% 
% % MLE of power law with cut off
% xmin = 30 ;
% sigIn = duration_phase ;
% sigIn = sigIn(sigIn>xmin) ;
% pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
% [lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
% L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;
% 
% c = histcounts(sigIn,linspace(1e1,1e5,1000),'Normalization','pdf');
% % loglog(x,c,'.','markersize',20)
% hold on
% loglog(x,pf_exp(x, lambdaHat3));
% 
% xlim([20 3e3])
% ylim([1e-5 1])
% legend('original data','exponential fitting')

% alpha stable fit

sigInProp = interval_phase ;
% pd = fitdist(sigInProp','exponential') ;

[n,x] = histcounts(sigInProp,400,'normalization','pdf') ;
semilogy((x(2:end)+x(1:end-1))/2,n,'.','MarkerSize',12)
hold on
xNew = linspace((x(2)+x(1))/2,(x(end)+x(end-1))/2,1000) ;
y = pdf(pd,xNew) ;
% semilogy(xNew,y,'r')
title('Interval distribution of randomized data')
xlabel('Time (ms)')
ylabel('probability')
xlim([20 4.5e3])   
ylim([1e-5 1e-1])
% legend('original data','\alpha stable distribution fit')
%% BIC


