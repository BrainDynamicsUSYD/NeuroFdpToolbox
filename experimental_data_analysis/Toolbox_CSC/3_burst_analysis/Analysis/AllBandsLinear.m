%% linear-spacing 
DurationFull = [] ;
load FullBurstYifanma027_032_Band1_3HzOctober24_15:12.mat Duration
% load FullBurstYifanma027_032_Band1_10HzOctober26_09:10.mat Duration	instantScale
instScale{1} = instantScale ;
DurationFull{1} = Duration ;

load FullBurstYifanma027_032_Band5_7HzOctober24_15:12.mat Duration instantScale
% load FullBurstYifanma027_032_Band11_20HzOctober26_09:13.mat Duration instantScale
instScale{2} = instantScale ;
DurationFull{2} = Duration ;

load FullBurstYifanma027_032_Band9_11HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band21_30HzOctober26_08:12.mat Duration instantScale
instScale{3} = instantScale ;
DurationFull{3} = Duration ;

load FullBurstYifanma027_032_Band13_15HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band31_40HzOctober26_08:12.mat	Duration instantScale
instScale{4} = instantScale ;
DurationFull{4} = Duration ;

load FullBurstYifanma027_032_Band17_19HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band41_50HzOctober26_08:12.mat	Duration instantScale
instScale{5} = instantScale ;
DurationFull{5} = Duration ;

load FullBurstYifanma027_032_Band21_23HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band51_60HzOctober26_08:13.mat	Duration instantScale
instScale{6} = instantScale ;
DurationFull{6} = Duration ;

load FullBurstYifanma027_032_Band25_27HzOctober24_14:12.mat Duration instantScale
% load FullBurstYifanma027_032_Band61_70HzOctober26_08:12.mat	Duration instantScale
instScale{7} = instantScale ;
DurationFull{7} = Duration ;

load FullBurstYifanma027_032_Band29_31HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band71_80HzOctober26_08:12.mat	Duration instantScale
instScale{8} = instantScale ;
DurationFull{8} = Duration ;

load FullBurstYifanma027_032_Band33_35HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band81_90HzOctober26_08:12.mat Duration instantScale
instScale{9} = instantScale ;
DurationFull{9} = Duration ;

load FullBurstYifanma027_032_Band37_39HzOctober24_14:13.mat Duration instantScale
% load FullBurstYifanma027_032_Band91_100HzOctober26_08:12.mat Duration instantScale
instScale{10} = instantScale ;
DurationFull{10} = Duration ;

load FullBurstYifanma027_032_Band41_43HzOctober24_15:10.mat Duration instantScale
DurationFull{11} = Duration ;
instScale{11} = instantScale ;

load FullBurstYifanma027_032_Band45_47HzOctober24_15:10.mat Duration instantScale
DurationFull{12} = Duration ;
instScale{12} = instantScale ;


%%
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'y-.'}] ;

for iBand = 1:10
    [n,x] = histcounts(DurationFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('2Hz','6Hz','10Hz','14Hz','18Hz','22Hz','26Hz',...
    '30Hz','34Hz','38Hz','42Hz','46Hz')
    
    
meanDu = [] ;
figure;
for iBand = 1:10
    meanDu(iBand) = mean(DurationFull{iBand}) ;
end
plot(meanDu,'o')

%% instant scale
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'},{'c-.'}] ;
meanIS = [] ;

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

figure;
plot(maxIS,'o')
semilogy(maxIS,'o')
ylabel('max burst size')
xlabel('frequency band index')
%% linear-spacing 
% load FullBurstYifanma027_032_Band1_3HzOctober24_15:12.mat patternScale
load FullBurstYifanma027_032_Band1_10HzOctober26_09:10.mat patternScale	
SizeFull{1} = patternScale ;

% load FullBurstYifanma027_032_Band5_7HzOctober24_15:12.mat patternScale
load FullBurstYifanma027_032_Band11_20HzOctober26_09:13.mat patternScale
SizeFull{2} = patternScale ;

% load FullBurstYifanma027_032_Band9_11HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band21_30HzOctober26_08:12.mat patternScale
SizeFull{3} = patternScale ;

% load FullBurstYifanma027_032_Band13_15HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band31_40HzOctober26_08:12.mat	patternScale
SizeFull{4} = patternScale ;

% load FullBurstYifanma027_032_Band17_19HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band41_50HzOctober26_08:12.mat	patternScale
SizeFull{5} = patternScale ;

% load FullBurstYifanma027_032_Band21_23HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band51_60HzOctober26_08:13.mat	patternScale
SizeFull{6} = patternScale ;

% load FullBurstYifanma027_032_Band25_27HzOctober24_14:12.mat patternScale
load FullBurstYifanma027_032_Band61_70HzOctober26_08:12.mat	patternScale
SizeFull{7} = patternScale ;

% load FullBurstYifanma027_032_Band29_31HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band71_80HzOctober26_08:12.mat	patternScale
SizeFull{8} = patternScale ;

% load FullBurstYifanma027_032_Band33_35HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band81_90HzOctober26_08:12.mat patternScale
SizeFull{9} = patternScale ;

% load FullBurstYifanma027_032_Band37_39HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band91_100HzOctober26_08:12.mat patternScale
SizeFull{10} = patternScale ;

% load FullBurstYifanma027_032_Band41_43HzOctober24_15:10.mat patternScale
% SizeFull{11} = patternScale ;

% load FullBurstYifanma027_032_Band45_47HzOctober24_15:10.mat patternScale
% SizeFull{12} = patternScale ;


%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'}] ;

for iBand = 1:10
    [n,x] = histcounts(SizeFull{iBand}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1-2Hz','3-4Hz','5-6Hz','7-8Hz','9-10Hz','11-13Hz','14-21Hz',...
    '22Hz-39Hz','30-40Hz','50-60Hz','70-80Hz')
    
meanSize = [] ;
figure;
for iBand = 1:10
    meanSize(iBand) = mean(SizeFull{iBand}) ;
end

plot(5:10:100,meanSize)


%% instant-spacing 
% load FullBurstYifanma027_032_Band1_3HzOctober24_15:12.mat patternScale
load FullBurstYifanma027_032_Band1_10HzOctober26_09:10.mat patternScale	
SizeFull{1} = patternScale ;

% load FullBurstYifanma027_032_Band5_7HzOctober24_15:12.mat patternScale
load FullBurstYifanma027_032_Band11_20HzOctober26_09:13.mat patternScale
SizeFull{2} = patternScale ;

% load FullBurstYifanma027_032_Band9_11HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band21_30HzOctober26_08:12.mat patternScale
SizeFull{3} = patternScale ;

% load FullBurstYifanma027_032_Band13_15HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band31_40HzOctober26_08:12.mat	patternScale
SizeFull{4} = patternScale ;

% load FullBurstYifanma027_032_Band17_19HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band41_50HzOctober26_08:12.mat	patternScale
SizeFull{5} = patternScale ;

% load FullBurstYifanma027_032_Band21_23HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band51_60HzOctober26_08:13.mat	patternScale
SizeFull{6} = patternScale ;

% load FullBurstYifanma027_032_Band25_27HzOctober24_14:12.mat patternScale
load FullBurstYifanma027_032_Band61_70HzOctober26_08:12.mat	patternScale
SizeFull{7} = patternScale ;

% load FullBurstYifanma027_032_Band29_31HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band71_80HzOctober26_08:12.mat	patternScale
SizeFull{8} = patternScale ;

% load FullBurstYifanma027_032_Band33_35HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band81_90HzOctober26_08:12.mat patternScale
SizeFull{9} = patternScale ;

% load FullBurstYifanma027_032_Band37_39HzOctober24_14:13.mat patternScale
load FullBurstYifanma027_032_Band91_100HzOctober26_08:12.mat patternScale
SizeFull{10} = patternScale ;

% load FullBurstYifanma027_032_Band41_43HzOctober24_15:10.mat patternScale
% SizeFull{11} = patternScale ;

% load FullBurstYifanma027_032_Band45_47HzOctober24_15:10.mat patternScale
% SizeFull{12} = patternScale ;