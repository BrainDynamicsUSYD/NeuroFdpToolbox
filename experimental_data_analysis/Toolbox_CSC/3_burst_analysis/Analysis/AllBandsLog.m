%% log-spacing 
load FullBurstYifanma027_032_Band127_129HzOctober24_10:31.mat Duration 
DurationFull{1} = Duration ;

load FullBurstYifanma027_032_Band63_65HzOctober24_10:31.mat Duration
DurationFull{2} = Duration ;

load FullBurstYifanma027_032_Band31_33HzOctober24_10:31.mat Duration
DurationFull{3} = Duration ;

load FullBurstYifanma027_032_Band15_17HzOctober24_10:31.mat Duration
DurationFull{4} = Duration ;

load FullBurstYifanma027_032_Band7_9HzOctober24_10:31.mat Duration
DurationFull{5} = Duration ;

load FullBurstYifanma027_032_Band3_5HzOctober24_11:31.mat Duration
DurationFull{6} = Duration ;

load FullBurstYifanma027_032_Band1_3HzOctober24_11:30.mat Duration
DurationFull{7} = Duration ;


%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'}] ;

for iBand = 1:7
    [n,x] = histcounts(DurationFull{7-iBand+1}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1-3Hz','3-5Hz','7-9Hz','15-17Hz','31-33Hz','63-65Hz',...
    '127Hz-129Hz')
    
    
meanDu = [] ;
figure;
for iBand = 1:7
    meanDu(iBand) = mean(DurationFull{7-iBand+1}) ;
end
plot(meanDu)



%% log-spacing 
load FullBurstYifanma027_032_Band127_129HzOctober24_10:31.mat patternScale 
SizeFull{1} = patternScale ;

load FullBurstYifanma027_032_Band63_65HzOctober24_10:31.mat patternScale
SizeFull{2} = patternScale ;

load FullBurstYifanma027_032_Band31_33HzOctober24_10:31.mat patternScale
SizeFull{3} = patternScale ;

load FullBurstYifanma027_032_Band15_17HzOctober24_10:31.mat patternScale
SizeFull{4} = patternScale ;

load FullBurstYifanma027_032_Band7_9HzOctober24_10:31.mat patternScale
SizeFull{5} = patternScale ;

load FullBurstYifanma027_032_Band3_5HzOctober24_11:31.mat patternScale
SizeFull{6} = patternScale ;

load FullBurstYifanma027_032_Band1_3HzOctober24_11:30.mat patternScale
SizeFull{7} = patternScale ;

%%
close all
figure;
shape = [{'r'},{'b'},{'k'},{'y'},{'g'},{'c'},{'m'},{'r-.'},{'k-.'},{'g-.'},{'b-.'}] ;

for iBand = 1:7
    [n,x] = histcounts(SizeFull{7-iBand+1}) ;
    loglog((x(1:end-1)+x(2:end))/2, n,shape{iBand}) ;
    hold on
end
legend('1-3Hz','3-5Hz','7-9Hz','15-17Hz','31-33Hz','63-65Hz',...
    '127Hz-129Hz')
    
meanSize = [] ;
figure;
for iBand = 1:7
    meanSize(iBand) = mean(SizeFull{7-iBand+1}) ;
end

plot(2.^(1:7),meanSize)

%% The location of bursts
load FullBurstYifanma027_032_Band1_3HzOctober24_11:30.mat rangeFrame patternScale instantScale
Delta_range = rangeFrame ;
Delta_size = patternScale ;
Delta_instantS = instantScale ;

load FullBurstYifanma027_032_Band3_5HzOctober24_11:31.mat rangeFrame patternScale instantScale
Theta_range = rangeFrame ;
Theta_size = patternScale ;
Theta_instantS = instantScale ;

load FullBurstYifanma027_032_Band7_9HzOctober24_10:31.mat rangeFrame patternScale instantScale
Alpha_range = rangeFrame ;
Alpha_size = patternScale ;
Alpha_instantS = instantScale ; 

load FullBurstYifanma027_032_Band15_17HzOctober24_10:31.mat rangeFrame patternScale instantScale
Beta_range = rangeFrame ;
Beta_size = patternScale ;
Beta_instantS = instantScale ; 

load FullBurstYifanma027_032_Band31_33HzOctober24_10:31.mat rangeFrame patternScale instantScale
Gamma1_range = rangeFrame ;
Gamma1_size = patternScale ;
Gamma1_instantS = instantScale ; 

load FullBurstYifanma027_032_Band63_65HzOctober24_10:31.mat rangeFrame patternScale instantScale
Gamma2_range = rangeFrame ;
Gamma2_size = patternScale ;
Gamma2_instantS = instantScale ; 

load FullBurstYifanma027_032_Band127_129HzOctober24_10:31.mat rangeFrame patternScale instantScale
Gamma3_range = rangeFrame ;
Gamma3_size = patternScale ;
Gamma3_instantS = instantScale ; 

%% Delta1 Gamma1
fsTemporal = 1.013e3 ;

Delta_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Delta_range,1)
    burstRange = Delta_range(iBurst,1):Delta_range(iBurst,2) ;
    Delta_burst(burstRange) = Delta_burst(burstRange)+Delta_instantS{iBurst} ;
end

Theta_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Theta_range,1)
    burstRange = Theta_range(iBurst,1):Theta_range(iBurst,2) ;
    Theta_burst(burstRange) = Theta_burst(burstRange)+Theta_instantS{iBurst} ;
end

Alpha_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Alpha_range,1)
    burstRange = Alpha_range(iBurst,1):Alpha_range(iBurst,2) ;
    Alpha_burst(burstRange) = Alpha_burst(burstRange)+Alpha_instantS{iBurst} ;
end

Beta_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Beta_range,1)
    burstRange = Beta_range(iBurst,1):Beta_range(iBurst,2) ;
    Beta_burst(burstRange) = Beta_burst(burstRange)+Beta_instantS{iBurst} ;
end

Gamma1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Gamma1_range,1)
    burstRange = Gamma1_range(iBurst,1):Gamma1_range(iBurst,2) ;
    Gamma1_burst(burstRange) = Gamma1_burst(burstRange)+Gamma1_instantS{iBurst} ;
end

figure; plot(Delta_burst); 
hold on;plot(Theta_burst-400,'r')
hold on;plot(Alpha_burst-800,'r')
hold on;plot(Beta_burst-1200,'r')
hold on;plot(Gamma1_burst-1600,'r')
xlabel('time (ms)')
ylabel('size')
title('temporal overalp of Delta burst and Gamma burst')
legend('Delta (1-2Hz)','Gamma (30-40Hz)')

