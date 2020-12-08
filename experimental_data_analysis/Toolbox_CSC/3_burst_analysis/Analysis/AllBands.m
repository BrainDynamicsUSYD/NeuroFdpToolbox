%%
load FullBurstYifanma027_032_Band0.5_80HzOctober22_11:56.mat Duration ...
    instantScale
DurationFull{1} = Duration ;
instScale{1} = instantScale ;

load FullBurstYifanma027_032_Band1_2HzOctober22_10:49.mat Duration ...
    instantScale
DurationFull{2} = Duration ;
instScale{2} = instantScale ;

load FullBurstYifanma027_032_Band3_4HzOctober22_10:50.mat Duration ...
    instantScale
DurationFull{3} = Duration ;
instScale{3} = instantScale ;

load FullBurstYifanma027_032_Band5_6HzOctober22_10:50.mat Duration ...
    instantScale
DurationFull{4} = Duration ;
instScale{4} = instantScale ;

load FullBurstYifanma027_032_Band7_8HzOctober22_10:51.mat Duration ...
    instantScale
DurationFull{5} = Duration ;
instScale{5} = instantScale ;

load FullBurstYifanma027_032_Band9_10HzOctober22_10:51.mat Duration ...
    instantScale
DurationFull{6} = Duration ;
instScale{6} = instantScale ;

load FullBurstYifanma027_032_Band11_13HzOctober22_11:52.mat Duration ...
    instantScale	  
DurationFull{7} = Duration ;
instScale{7} = instantScale ;

load FullBurstYifanma027_032_Band14_21HzOctober22_11:55.mat Duration ...
    instantScale	 
DurationFull{8} = Duration ;
instScale{8} = instantScale ;

load FullBurstYifanma027_032_Band22_29HzOctober22_11:52.mat Duration ...
    instantScale	 
DurationFull{9} = Duration ;
instScale{9} = instantScale ;

load FullBurstYifanma027_032_Band30_40HzOctober22_10:57.mat Duration ...
    instantScale	 
DurationFull{10} = Duration ;
instScale{10} = instantScale ;

load FullBurstYifanma027_032_Band50_60HzOctober22_10:56.mat Duration ...
    instantScale
DurationFull{11} = Duration ;
instScale{11} = instantScale ;

load FullBurstYifanma027_032_Band70_80HzOctober22_10:56.mat Duration ...
    instantScale
DurationFull{12} = Duration ;
instScale{12} = instantScale ;

%%
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'}] ;

for iBand = 1:11
    [n,x] = histcounts(DurationFull{iBand+1}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1-2Hz','3-4Hz','5-6Hz','7-8Hz','9-10Hz','11-13Hz','14-21Hz',...
    '22Hz-39Hz','30-40Hz','50-60Hz','70-80Hz')
    
figure;
for iBand = 1:11
    meanIS(iBand) = mean(DurationFull{iBand+1}) ;
end
plot(meanIS)

%% instant scale
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;
meanIS = [] ;

for iBand = 1:12
    tempScale = instScale{iBand} ;
    instS = [] ;
    for iBurst = 1:length(tempScale)
        instS = [instS;tempScale{iBurst}] ;        
    end  
    meanIS(iBand) = mean(instS) ;
    [n,x] = histcounts(instS) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
figure;
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst size')
xlabel('frequency band index')

%%
load FullBurstYifanma027_032_Band0.5_80HzOctober22_11:56.mat patternScale
SizeFull{1} = patternScale ;

load FullBurstYifanma027_032_Band1_2HzOctober22_10:49.mat patternScale	 
SizeFull{2} = patternScale ;

load FullBurstYifanma027_032_Band3_4HzOctober22_10:50.mat patternScale
SizeFull{3} = patternScale ;

load FullBurstYifanma027_032_Band5_6HzOctober22_10:50.mat patternScale
SizeFull{4} = patternScale ;

load FullBurstYifanma027_032_Band7_8HzOctober22_10:51.mat patternScale
SizeFull{5} = patternScale ;

load FullBurstYifanma027_032_Band9_10HzOctober22_10:51.mat patternScale
SizeFull{6} = patternScale ;

load FullBurstYifanma027_032_Band11_13HzOctober22_11:52.mat patternScale	  
SizeFull{7} = patternScale ;

load FullBurstYifanma027_032_Band14_21HzOctober22_11:55.mat patternScale	 
SizeFull{8} = patternScale ;

load FullBurstYifanma027_032_Band22_29HzOctober22_11:52.mat patternScale	 
SizeFull{9} = patternScale ;

load FullBurstYifanma027_032_Band30_40HzOctober22_10:57.mat patternScale	 
SizeFull{10} = patternScale ;

load FullBurstYifanma027_032_Band50_60HzOctober22_10:56.mat patternScale
SizeFull{11} = patternScale ;

load FullBurstYifanma027_032_Band70_80HzOctober22_10:56.mat patternScale
SizeFull{12} = patternScale ;


%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'}] ;

for iBand = 1:11
    [n,x] = histcounts(SizeFull{iBand+1}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1-2Hz','3-4Hz','5-6Hz','7-8Hz','9-10Hz','11-13Hz','14-21Hz',...
    '22Hz-39Hz','30-40Hz','50-60Hz','70-80Hz')
    

figure;
for iBand = 1:11
    meanSize(iBand) = mean(SizeFull{iBand+1}) ;
end
plot(meanSize)



%% according to neuroscience
DurationFull = [] ;
load FullBurstYifanma027_032_Band0.5_1.5HzJanuary10_09:21.mat Duration ...
    patternScale instantScale
DurationFull{1} = Duration ;
SizeFull{1} = patternScale ;
instScale{1} = instantScale ;

load FullBurstYifanma027_032_Band1.5_2.5HzJanuary10_09:20.mat Duration ...
    patternScale instantScale
DurationFull{2} = Duration ;
SizeFull{2} = patternScale ;
instScale{2} = instantScale ;

load FullBurstYifanma027_032_Band2.5_4HzJanuary10_08:21.mat Duration ...
    instantScale patternScale	
DurationFull{3} = Duration ;
SizeFull{3} = patternScale ;
instScale{3} = instantScale ;

load FullBurstYifanma027_032_Band4_6HzJanuary10_08:21.mat Duration ...
    instantScale patternScale
DurationFull{4} = Duration ;
SizeFull{4} = patternScale ;
instScale{4} = instantScale ;

load FullBurstYifanma027_032_Band6_8HzJanuary10_08:21.mat Duration ...
    instantScale patternScale
DurationFull{5} = Duration ;
SizeFull{5} = patternScale ;
instScale{5} = instantScale ;

load FullBurstYifanma027_032_Band8_10HzJanuary10_08:21.mat Duration ...
    instantScale patternScale
DurationFull{6} = Duration ;
SizeFull{6} = patternScale ;
instScale{6} = instantScale ;

load FullBurstYifanma027_032_Band10_13HzJanuary10_08:21.mat	Duration ...
    instantScale patternScale
DurationFull{7} = Duration ;
SizeFull{7} = patternScale ;
instScale{7} = instantScale ;

load FullBurstYifanma027_032_Band13_20HzJanuary10_08:22.mat	Duration ...
    instantScale patternScale
DurationFull{8} = Duration ;
SizeFull{8} = patternScale ;
instScale{8} = instantScale ;

load FullBurstYifanma027_032_Band20_30HzJanuary10_08:23.mat Duration ...
    instantScale patternScale
DurationFull{9} = Duration ;
SizeFull{9} = patternScale ;
instScale{9} = instantScale ;

load FullBurstYifanma027_032_Band30_50HzJanuary10_08:25.mat Duration ...
    instantScale patternScale
DurationFull{10} = Duration ;
SizeFull{10} = patternScale ;
instScale{10} = instantScale ;

load FullBurstYifanma027_032_Band50_80HzJanuary10_09:27.mat Duration patternScale
DurationFull{11} = Duration ;
SizeFull{11} = patternScale ;
instScale{11} = instantScale ;

load FullBurstYifanma027_032_Band80_120HzJanuary10_09:29.mat Duration patternScale
DurationFull{12} = Duration ;
SizeFull{12} = patternScale ;
instScale{12} = instantScale ;

%%
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:12
    [n,x] = histcounts(DurationFull{iBand},'normalization','pdf') ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}, 'markersize',12) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:12
    meanIS(iBand) = mean(DurationFull{iBand}) ;
end
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst duration')
xlabel('frequency band index')
%%
figure;
titleName =  [{'Delta1'},{'Delta2'},{'Delta3'},{'Theta1'},{'Theta2'},{'Alpha1'},{'Alpha2'},...
    {'Beta1'},{'Beta2'},{'Gamma1'},{'Gamma2'},{'Gamma3'}] ;
for iBand = 1:12
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(DurationFull{iBand},10+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['duration dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e1, 4e3])
    xlabel('time')
    ylabel('probability')
end
% legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
%    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')

%% size of neuroscience spacing
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:12
    [n,x] = histcounts(SizeFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:12
    meanIS(iBand) = mean(SizeFull{iBand}) ;
end
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst size')
xlabel('frequency band index')
title('mean size of amp. pattern across scales')

%
figure;
titleName =  [{'Delta1'},{'Delta2'},{'Delta3'},{'Theta1'},{'Theta2'},{'Alpha1'},{'Alpha2'},...
    {'Beta1'},{'Beta2'},{'Gamma1'},{'Gamma2'},{'Gamma3'}] ;
for iBand = 1:12
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(SizeFull{iBand},10+iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e2, 8e5])
    xlabel('size')
    ylabel('probability')
end


% legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
%    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')

%% instant scale
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;
meanIS = [] ;

for iBand = 1:12
    tempScale = instScale{iBand} ;
    instS = [] ;
    for iBurst = 1:length(tempScale)
        instS = [instS;tempScale{iBurst}] ;        
    end  
    meanIS(iBand) = mean(instS) ;
    [n,x] = histcounts(instS) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
figure;
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst size')
xlabel('frequency band index')


for iBand = 1:12
    tempScale = instScale{iBand} ;
    instS = [] ;
    maxIStemp = 0 ;
    for iBurst = 1:length(tempScale)
        instS = [instS;tempScale{iBurst}] ; 
        maxIStemp = maxIStemp + max(tempScale{iBurst}) ;
    end  
    meanIS(iBand) = mean(instS) ;
    maxIS(iBand) = maxIStemp/iBurst ;
%     [n,x] = histcounts(instS) ;
%     loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
%     hold on
    
    
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(SizeFull{iBand},10+iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e2, 8e5])
    xlabel('size')
    ylabel('probability')
    
end


figure;
plot(maxIS,'o')
semilogy(maxIS,'o')
ylabel('max burst size')
xlabel('frequency band index')

%
figure;
titleName =  [{'Delta1'},{'Delta2'},{'Delta3'},{'Theta1'},{'Theta2'},{'Alpha1'},{'Alpha2'},...
    {'Beta1'},{'Beta2'},{'Gamma1'},{'Gamma2'},{'Gamma3'}] ;
for iBand = 1:12
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(SizeFull{iBand},10+iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e2, 8e5])
    xlabel('size')
    ylabel('probability')
end



%% according to logspace (log bandwidth)
DurationFull = [] ;
SizeFull = [] ;

load FullBurstYifanma027_032_Band0.75_1.25HzJanuary10_09:31.mat Duration patternScale
DurationFull{1} = Duration ;
SizeFull{1} = patternScale ;

load FullBurstYifanma027_032_Band1.5_2.5HzJanuary10_09:31.mat Duration patternScale  
DurationFull{2} = Duration ;
SizeFull{2} = patternScale ;

load FullBurstYifanma027_032_Band3_5HzJanuary10_08:32.mat Duration patternScale
DurationFull{3} = Duration ;
SizeFull{3} = patternScale ;

load FullBurstYifanma027_032_Band6_10HzJanuary10_08:32.mat Duration patternScale
DurationFull{4} = Duration ;
SizeFull{4} = patternScale ;

load FullBurstYifanma027_032_Band12_20HzJanuary10_08:33.mat Duration patternScale
DurationFull{5} = Duration ;
SizeFull{5} = patternScale ;

load FullBurstYifanma027_032_Band24_40HzJanuary10_08:36.mat Duration patternScale
DurationFull{6} = Duration ;
SizeFull{6} = patternScale ;

load FullBurstYifanma027_032_Band48_80HzJanuary10_08:38.mat Duration patternScale
DurationFull{7} = Duration ;
SizeFull{7} = patternScale ;

load FullBurstYifanma027_032_Band96_160HzJanuary10_08:41.mat Duration patternScale
DurationFull{8} = Duration ;
SizeFull{8} = patternScale ;

%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:8
    [n,x] = histcounts(DurationFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1Hz','2Hz','4Hz','8Hz','16Hz','32Hz','64Hz',...
    '128Hz')
    
    
meanIS = [] ;
figure;
for iBand = 1:8
    meanIS(iBand) = mean(DurationFull{iBand}) ;
end
% plot(meanDu,'o')
loglog(2.^(0:7),meanIS(1:end),'o')
ylabel('mean burst duration')
xlabel('frequency')

%
figure;
titleName =  [{'1Hz'},{'2Hz'},{'4Hz'},{'8Hz'},{'16Hz'},{'32Hz'},{'64Hz'},...
    {'128Hz'}] ;
for iBand = 1:8
    subplot(8,1,(iBand))
     histogram(DurationFull{iBand},15+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['duration dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e1, 1e4])
    xlabel('time (ms)')
    ylabel('probability')
end

%% size log-spacing
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:8
    [n,x] = histcounts(SizeFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:8
    meanIS(iBand) = mean(SizeFull{iBand}) ;
end
plot(meanIS,'o')
loglog(2.^(0:7),meanIS,'o')
ylabel('mean burst size')
xlabel('frequency (Hz)')
title('mean size of amp. pattern across scales')

%
figure;
titleName =  [{'1Hz'},{'2Hz'},{'4Hz'},{'8Hz'},{'16Hz'},{'32Hz'},{'64Hz'},...
    {'128Hz'}] ;
for iBand = 1:8
    subplot(8,1,(iBand))
     histogram(SizeFull{iBand},10+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e0, 1.2e6])
    xlabel('size')
    ylabel('probability')
end


%% test data (swap sequence)
DurationFull = [] ;

load FullBurstYifanma027_032_Band1_10HzJanuary10_10:47.mat Duration	
DurationFull{1} = Duration ;

load FullBurstYifanma027_032_Band11_18HzJanuary10_10:47.mat	Duration
DurationFull{2} = Duration ;

load FullBurstYifanma027_032_Band19_25HzJanuary10_09:47.mat	Duration
DurationFull{3} = Duration ;

load FullBurstYifanma027_032_Band26_30HzJanuary10_09:47.mat	Duration
DurationFull{4} = Duration ;

load FullBurstYifanma027_032_Band31_34HzJanuary10_09:47.mat Duration
DurationFull{5} = Duration ;

load FullBurstYifanma027_032_Band35_37HzJanuary10_09:47.mat Duration
DurationFull{6} = Duration ;

load FullBurstYifanma027_032_Band38_39HzJanuary10_09:46.mat Duration
DurationFull{7} = Duration ;

load FullBurstYifanma027_032_Band40_41HzJanuary10_09:46.mat Duration
DurationFull{8} = Duration ;

%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:8
    [n,x] = histcounts(DurationFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1Hz','2Hz','4Hz','8Hz','16Hz','32Hz','64Hz',...
    '128Hz')
    
    
meanIS = [] ;
figure;
for iBand = 1:8
    meanIS(iBand) = mean(DurationFull{iBand}) ;
end
% plot(meanDu,'o')
loglog(2.^(0:7),meanIS(1:end),'o')
ylabel('mean burst duration')
xlabel('frequency')

%% surrogate data
DurationFull = [] ;
SizeFull = [] ;

load FullBurstYifanma027_032_Band0.5_1.5HzJanuary11_10:54.mat patternScale
DurationFull{1} = Duration ;
SizeFull{1} = patternScale ;

load FullBurstYifanma027_032_Band1.5_2.5HzJanuary11_10:54.mat patternScale
DurationFull{2} = Duration ;
SizeFull{2} = patternScale ;

load FullBurstYifanma027_032_Band2.5_4HzJanuary11_09:54.mat patternScale	  
DurationFull{3} = Duration ;
SizeFull{3} = patternScale ;

load FullBurstYifanma027_032_Band4_6HzJanuary11_09:54.mat patternScale
DurationFull{4} = Duration ;
SizeFull{4} = patternScale ;

load FullBurstYifanma027_032_Band6_8HzJanuary11_09:54.mat patternScale
DurationFull{5} = Duration ;
SizeFull{5} = patternScale ;

load FullBurstYifanma027_032_Band8_10HzJanuary11_09:54.mat patternScale
DurationFull{6} = Duration ;
SizeFull{6} = patternScale ;

load FullBurstYifanma027_032_Band10_13HzJanuary11_09:54.mat	patternScale  
DurationFull{7} = Duration ;
SizeFull{7} = patternScale ;

load FullBurstYifanma027_032_Band13_20HzJanuary11_09:55.mat	patternScale
DurationFull{8} = Duration ;
SizeFull{8} = patternScale ;

load FullBurstYifanma027_032_Band20_30HzJanuary11_09:55.mat	patternScale 
DurationFull{9} = Duration ;
SizeFull{9} = patternScale ;

load FullBurstYifanma027_032_Band30_50HzJanuary11_09:56.mat patternScale
DurationFull{10} = Duration ;
SizeFull{10} = patternScale ;

load FullBurstYifanma027_032_Band50_80HzJanuary11_10:54.mat patternScale
DurationFull{11} = Duration ;
SizeFull{11} = patternScale ;

load FullBurstYifanma027_032_Band80_120HzJanuary11_10:53.mat patternScale
DurationFull{12} = Duration ;
SizeFull{12} = patternScale ;


%%
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:12
    [n,x] = histcounts(DurationFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:12
    meanIS(iBand) = mean(DurationFull{iBand}) ;
end
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst duration')
xlabel('frequency band index')

%%
figure;
titleName =  [{'Delta1'},{'Delta2'},{'Delta3'},{'Theta1'},{'Theta2'},{'Alpha1'},{'Alpha2'},...
    {'Beta1'},{'Beta2'},{'Gamma1'},{'Gamma2'},{'Gamma3'}] ;
for iBand = 1:12
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(DurationFull{iBand},10+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['duration dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e1, 4e3])
    xlabel('time')
    ylabel('probability')
end
% legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
%    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')

%% size of neuroscience spacing
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:12
    [n,x] = histcounts(SizeFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:12
    meanIS(iBand) = mean(SizeFull{iBand}) ;
end
plot(meanIS,'o')
semilogy(meanIS,'o')
ylabel('mean burst size')
xlabel('frequency band index')
title('mean size of amp. pattern across scales')

%
figure;
titleName =  [{'Delta1'},{'Delta2'},{'Delta3'},{'Theta1'},{'Theta2'},{'Alpha1'},{'Alpha2'},...
    {'Beta1'},{'Beta2'},{'Gamma1'},{'Gamma2'},{'Gamma3'}] ;
for iBand = 1:12
    plotNum = [1:2:11,2:2:12] ;
    subplot(6,2,plotNum(iBand))
     histogram(SizeFull{iBand},10+iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e2, 2e5])
    xlabel('size')
    ylabel('probability')
end


% legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
%    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')


%% according to surrogate logspace (log bandwidth)
DurationFull = [] ;
SizeFull = [] ;

load FullBurstYifanma027_032_Band0.75_1.25HzJanuary14_15:45.mat Duration patternScale
DurationFull{1} = Duration ;
SizeFull{1} = patternScale ;

load FullBurstYifanma027_032_Band1.5_2.5HzJanuary14_15:45.mat Duration patternScale  
DurationFull{2} = Duration ;
SizeFull{2} = patternScale ;

load FullBurstYifanma027_032_Band3_5HzJanuary14_14:44.mat Duration patternScale
DurationFull{3} = Duration ;
SizeFull{3} = patternScale ;

load FullBurstYifanma027_032_Band6_10HzJanuary14_14:45.mat Duration patternScale
DurationFull{4} = Duration ;
SizeFull{4} = patternScale ;

load FullBurstYifanma027_032_Band12_20HzJanuary14_14:45.mat Duration patternScale
DurationFull{5} = Duration ;
SizeFull{5} = patternScale ;

load FullBurstYifanma027_032_Band24_40HzJanuary14_14:46.mat Duration patternScale
DurationFull{6} = Duration ;
SizeFull{6} = patternScale ;

load FullBurstYifanma027_032_Band48_80HzJanuary14_14:47.mat Duration patternScale
DurationFull{7} = Duration ;
SizeFull{7} = patternScale ;

load FullBurstYifanma027_032_Band96_160HzJanuary14_14:45.mat Duration patternScale
DurationFull{8} = Duration ;
SizeFull{8} = patternScale ;

%% surrogate log
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:8
    [n,x] = histcounts(DurationFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1Hz','2Hz','4Hz','8Hz','16Hz','32Hz','64Hz',...
    '128Hz')
    
    
meanIS = [] ;
figure;
for iBand = 1:8
    meanIS(iBand) = mean(DurationFull{iBand}) ;
end
% plot(meanDu,'o')
loglog(2.^(0:7),meanIS(1:end),'o')
ylabel('mean burst duration')
xlabel('frequency')

%
figure;
titleName =  [{'1Hz'},{'2Hz'},{'4Hz'},{'8Hz'},{'16Hz'},{'32Hz'},{'64Hz'},...
    {'128Hz'}] ;
for iBand = 1:8
    subplot(8,1,(iBand))
     histogram(DurationFull{iBand},15+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['duration dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e1, 1e4])
    xlabel('time (ms)')
    ylabel('probability')
end

%% size surrogate log-spacing
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;

for iBand = 1:8
    [n,x] = histcounts(SizeFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('Delta1','Delta2','Delta3','Theta1','Theta2','Alpha1','Alpha2',...
    'Beta1','Beta2','Gamma1','Gamma2','Gamma3')
    
    
meanIS = [] ;
figure;
for iBand = 1:8
    meanIS(iBand) = mean(SizeFull{iBand}) ;
end
plot(meanIS,'o')
loglog(2.^(0:7),meanIS,'o')
ylabel('mean burst size')
xlabel('frequency (Hz)')
title('mean size of amp. pattern across scales')

%
figure;
titleName =  [{'1Hz'},{'2Hz'},{'4Hz'},{'8Hz'},{'16Hz'},{'32Hz'},{'64Hz'},...
    {'128Hz'}] ;
for iBand = 1:8
    subplot(8,1,(iBand))
     histogram(SizeFull{iBand},10+2*iBand ,'normalization','pdf') ;
    % loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    title(['size dist. of ', titleName{iBand}, ' amp. pattern'])
    xlim([1e0, 1.2e6])
    xlabel('size')
    ylabel('probability')
end