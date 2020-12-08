%% study the overlap and pattern correlation between Gamma and Delta
%
load FullBurstYifanma027_032_Band1_2HzOctober11_19:31.mat CC burstIdx
Delta1_CC = CC ;
Delta1_burstIdx = burstIdx ;
load FullBurstYifanma027_032_Band3_4HzOctober11_19:31.mat CC burstIdx
Delta2_CC = CC ;
Delta2_burstIdx = burstIdx ;
load FullBurstYifanma027_032_Band30_40HzOctober11_19:38.mat CC burstIdx
Gamma1_CC = CC ;
Gamma1_burstIdx = burstIdx ;
load FullBurstYifanma027_032_Band50_60HzOctober11_19:37.mat CC burstIdx
Gamma2_CC = CC ;
Gamma2_burstIdx = burstIdx ;
load FullBurstYifanma027_032_Band70_80HzOctober11_19:38.mat CC burstIdx
Gamma3_CC = CC ;
Gamma3_burstIdx = burstIdx ;

%% Delta1 vs Gamma1
fsTemporal = 1.013e3 ;
DeltaBurst = Delta1_CC.PixelIdxList ;
GammaBurst = Gamma1_CC.PixelIdxList ;

activity = zeros(20,20,310*fsTemporal) ;

for iBurst = 1:length(DeltaBurst)
    activity(DeltaBurst{iBurst}) = activity(DeltaBurst{iBurst})+1 ;
end

for iBurst = 1:length(GammaBurst)
    activity(GammaBurst{iBurst}) = activity(GammaBurst{iBurst})+2 ;
end
%
no_overlap = length([find(activity==1);find(activity==2)]) ;
overlap = length(find(activity==3)) ;
Delta_space = length(find(activity==1)) ;
Gamma_space = length(find(activity==2)) ;

%% Delta2 vs Gamma1
fsTemporal = 1.013e3 ;
DeltaBurst = Delta2_CC.PixelIdxList ;
GammaBurst = Gamma1_CC.PixelIdxList ;

activity = zeros(20,20,310*fsTemporal) ;

for iBurst = 1:length(DeltaBurst)
    activity(DeltaBurst{iBurst}) = activity(DeltaBurst{iBurst})+1 ;
end

for iBurst = 1:length(GammaBurst)
    activity(GammaBurst{iBurst}) = activity(GammaBurst{iBurst})+1 ;
end

no_overlap2 = length(find(activity==1)) ;
overlap2 = length(find(activity==2)) ;

%% Delta1 vs Gamma2
fsTemporal = 1.013e3 ;
DeltaBurst = Delta1_CC.PixelIdxList ;
GammaBurst = Gamma2_CC.PixelIdxList ;

activity = zeros(20,20,310*fsTemporal) ;

for iBurst = 1:length(DeltaBurst)
    activity(DeltaBurst{iBurst}) = activity(DeltaBurst{iBurst})+1 ;
end

for iBurst = 1:length(GammaBurst)
    activity(GammaBurst{iBurst}) = activity(GammaBurst{iBurst})+1 ;
end

no_overlap3 = length(find(activity==1)) ;
overlap3 = length(find(activity==2)) ;

%% Delta1 vs Gamma3
fsTemporal = 1.013e3 ;
DeltaBurst = Delta1_CC.PixelIdxList ;
GammaBurst = Gamma3_CC.PixelIdxList ;

activity = zeros(20,20,310*fsTemporal) ;

for iBurst = 1:length(DeltaBurst)
    activity(DeltaBurst{iBurst}) = activity(DeltaBurst{iBurst})+1 ;
end

for iBurst = 1:length(GammaBurst)
    activity(GammaBurst{iBurst}) = activity(GammaBurst{iBurst})+1 ;
end

no_overlap4 = length(find(activity==1)) ;
overlap4 = length(find(activity==2)) ;

%% visualisation

vidTitle = 'Delta_Gamma_overlap3'
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = 20 ;
    open(vidObj);
    fig=figure ;
    set(gcf,'Position',[260 23 1159 926])
	timeCount = 0 ;
    
for iTime = 10200:20200
    imagesc(activity(:,:,iTime))
    caxis([0 3])
    colorbar
    title(['Delta:blue; Gamma: dark yellow; Overlap: yellow; time: ',num2str(iTime),'ms'])
    writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    cla
end
close(vidObj);

%% The location of bursts
load FullBurstYifanma027_032_Band1_2HzOctober11_19:31.mat rangeFrame patternScale instantScale
Delta1_range = rangeFrame ;
Delta1_size = patternScale ;
Delta1_instantS = instantScale ; 

load FullBurstYifanma027_032_Band3_4HzOctober11_19:31.mat rangeFrame patternScale instantScale
Delta2_range = rangeFrame ;
Delta2_size = patternScale ;
Delta2_instantS = instantScale ; 

load FullBurstYifanma027_032_Band30_40HzOctober11_19:38.mat rangeFrame patternScale instantScale
Gamma1_range = rangeFrame ;
Gamma1_size = patternScale ;
Gamma1_instantS = instantScale ; 

load FullBurstYifanma027_032_Band50_60HzOctober11_19:37.mat rangeFrame patternScale instantScale
Gamma2_range = rangeFrame ;
Gamma2_size = patternScale ;
Gamma2_instantS = instantScale ; 

load FullBurstYifanma027_032_Band70_80HzOctober11_19:38.mat rangeFrame patternScale instantScale
Gamma3_range = rangeFrame ;
Gamma3_size = patternScale ;
Gamma3_instantS = instantScale ; 


%% Delta1 Gamma1
Delta1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Delta1_range,1)
    burstRange = Delta1_range(iBurst,1):Delta1_range(iBurst,2) ;
    Delta1_burst(burstRange) = Delta1_burst(burstRange)+Delta1_instantS{iBurst} ;
end

Gamma1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Gamma1_range,1)
    burstRange = Gamma1_range(iBurst,1):Gamma1_range(iBurst,2) ;
    Gamma1_burst(burstRange) = Gamma1_burst(burstRange)+Gamma1_instantS{iBurst} ;
end

figure; plot(Delta1_burst); hold on;plot(Gamma1_burst-400,'r')
xlabel('time (ms)')
ylabel('size')
title('temporal overalp of Delta burst and Gamma burst')
legend('Delta (1-2Hz)','Gamma (30-40Hz)')

%% correlation
xcov(Delta1_burst-mean(Delta1_burst),Gamma1_burst-mean(Gamma1_burst)) ;
xcorr(Delta1_burst-mean(Delta1_burst),Gamma1_burst-mean(Gamma1_burst)) ;
figure;
plot(xcorr)

corr2(Delta1_burst,Gamma1_burst)

%% Delta2 Gamma1
Delta2_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Delta2_range,1)
    burstRange = Delta2_range(iBurst,1):Delta2_range(iBurst,2) ;
    Delta2_burst(burstRange) = Delta2_burst(burstRange)+Delta2_instantS{iBurst} ;
end

Gamma1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Gamma1_range,1)
    burstRange = Gamma1_range(iBurst,1):Gamma1_range(iBurst,2) ;
    Gamma1_burst(burstRange) = Gamma1_burst(burstRange)+Gamma1_instantS{iBurst} ;
end

figure; plot(Delta2_burst); hold on;plot(Gamma1_burst-400,'r')
xlabel('time (ms)')
ylabel('size')
title('temporal overalp of Delta burst and Gamma burst')
legend('Delta (3-4Hz)','Gamma (30-40Hz)')

%% Delta1 Gamma2
Delta1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Delta1_range,1)
    burstRange = Delta1_range(iBurst,1):Delta1_range(iBurst,2) ;
    Delta1_burst(burstRange) = Delta1_burst(burstRange)+Delta1_instantS{iBurst} ;
end

Gamma2_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Gamma2_range,1)
    burstRange = Gamma2_range(iBurst,1):Gamma2_range(iBurst,2) ;
    Gamma2_burst(burstRange) = Gamma2_burst(burstRange)+Gamma2_instantS{iBurst} ;
end

figure; plot(Delta1_burst); hold on;plot(Gamma2_burst-400,'r')
xlabel('time (ms)')
ylabel('size')
title('temporal overalp of Delta burst and Gamma burst')
legend('Delta (1-2Hz)','Gamma (50-60Hz)')

%% Delta1 Gamma3
Delta1_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Delta1_range,1)
    burstRange = Delta1_range(iBurst,1):Delta1_range(iBurst,2) ;
    Delta1_burst(burstRange) = Delta1_burst(burstRange)+Delta1_instantS{iBurst} ;
end

Gamma3_burst = zeros(310*fsTemporal,1) ;
for iBurst = 1:size(Gamma3_range,1)
    burstRange = Gamma3_range(iBurst,1):Gamma3_range(iBurst,2) ;
    Gamma3_burst(burstRange) = Gamma3_burst(burstRange)+Gamma3_instantS{iBurst} ;
end

figure; plot(Delta1_burst); hold on;plot(Gamma3_burst-400,'r')
xlabel('time (ms)')
ylabel('size')
title('temporal overalp of Delta burst and Gamma burst')
legend('Delta (1-2Hz)','Gamma (70-80)')

%% look at the speed of burst pattern




