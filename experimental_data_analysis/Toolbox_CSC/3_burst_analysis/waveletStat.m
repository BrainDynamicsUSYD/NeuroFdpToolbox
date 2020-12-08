%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation
clear ; 
close ;
clc ;
cd ..
%% add path
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/uimage']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt/NeuroPattToolbox-master']))
addpath(genpath([pwd,'/ToolOthers/UTILITIES_Paul']))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in data
dataFileName = 'ma027_032' ;
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing (optional)
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

% %% Find subband signals
% % apply bandpass filtering for Gamma signals
% subBand = [30,80] ;
% bandpassSig = find_bandpassSig(sigOri,subBand, fsTemporal,3) ;
% % Hilbert transform for analytic signals
% hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
% % find amplitdue of the analytic Gamma as the input
% sigIn = abs(squeeze((hilbertSig))) ;
Dfre = [] ;
Duration = [] ;
interval = [] ;
for iChannel = setdiff(1:100,badChannels)
    
    
    
sigIn = squeeze(sigOri(mod(iChannel-1,10)+1,ceil((iChannel-1)/10)+1,:)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% wName = 'cmor1.5-1' ;
% freqRange = [30, 80] ;
% scaleStep = 0.5 ;
% fc = centfrq(wName) ;
% scalerange = fc./(freqRange/fsTemporal) ;
% scales = scalerange(end):scaleStep:scalerange(1) ;            % 
% pseudoFreq = scal2frq(scales, wName, 1/fsTemporal) ;
% 
% wt = cwt( sigIn ,scales, wName  ) ;
% f1 = pseudoFreq ;

[wt2,f22,coi] = cwt(sigIn,fsTemporal,'VoicesPerOctave',30);
    for j = 1:length(coi)
        ind = find(f22<=coi(j));
        wt2(ind,j) = NaN;
    end
 c

%%

% CData = abs(wt);
% % CData = abs(wt2);
%     Y = prctile(CData(:),97); % 95
%     CData(CData < Y) = 0;
%     GreyImage1 = CData;
%     CData(CData>0) = 1;
%     binaryImage1 = CData;
    
    CData2 = abs(wt2);
% CData = abs(wt2);
    Y = prctile(CData2(:),95); % 95
    CData2(CData2 < Y) = 0;
    GreyImage2 = CData2;
    CData2(CData2>0) = 1;
    binaryImage2 = CData2;
    
    %%
% close all
% timePeriod = 230000:234000 ;
% figure;imagesc(timePeriod,f1,abs(wt(:,timePeriod)))
% set(gca,'YDir','normal')
% figure;imagesc(timePeriod,f2,abs(wt2(:,timePeriod)))
% set(gca,'YDir','normal')
% figure;imagesc(timePeriod,f1,binaryImage1(:,timePeriod))
% set(gca,'YDir','normal')
% figure;imagesc(timePeriod,f2,binaryImage2(:,timePeriod))
% set(gca,'YDir','normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find subband signals
% % apply bandpass filtering for Gamma signals
% subBand = [30,80] ;
% bandpassSig = find_bandpassSig(sigIn',subBand, fsTemporal,2,[],[],8,0) ;
% % Hilbert transform for analytic signals
% hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4,[],0) ;
% % find amplitdue of the analytic Gamma as the input
% sigHilbert = abs(squeeze((hilbertSig))) ;

%% Find Gamma burst in the time domain for fig.1
% stdVal = 3;
% GammaBurstEvent = find_Burst_1D(sigHilbert',fsTemporal,0,[],...
%     1,stdVal) ;
%%
% figure;
% plot(timePeriod,squeeze(bandpassSig(1,1,1,timePeriod)))
% hold on
% plot(timePeriod,sigHilbert(timePeriod))
% hold on
% plot(timePeriod,GammaBurstEvent.is_burst(timePeriod)*max(sigHilbert(:)))
% xlim([timePeriod(1) timePeriod(end)])

%%
freq = f2 ;
imageIn = binaryImage2 ;
GreyImage = GreyImage2 ;
tic
    CC = bwconncomp(imageIn) ;
     B = regionprops(CC,'BoundingBox'); % Smallest rectangle containing the region
    Boundary = cat(1, B.BoundingBox); % [ul_corner width], ul_corner specifies
    % the upper-left corner [x y z ...]. width specifies the width along each dimension [x_width y_width ...]
    stTime =  Boundary(:,1)/fsTemporal*1000; % milisecond
    duTime = Boundary(:,3)/fsTemporal*1000 ;
    endTime = stTime + duTime;    
    cFreq = [];
    iRegion = 2;
    j = 0;
    while iRegion < length(stTime)
        if stTime(iRegion+1) > endTime(iRegion)  %(size(CC.PixelIdxList{iRegion},1)>200)
%              certain = GreyImage(CC.PixelIdxList{iRegion});
%              [I,~] = ind2sub(size(GreyImage),CC.PixelIdxList{iRegion}(find(certain == max(certain))));
           instantBinary = imageIn(:,Boundary(iRegion,1):Boundary(iRegion+j,1)+Boundary(iRegion+j,3));
           instantPattern = GreyImage(:,Boundary(iRegion,1):Boundary(iRegion+j,1)+Boundary(iRegion+j,3)); 
        else
            while iRegion+1+j <= length(stTime) && stTime(iRegion+1+j) <= endTime(iRegion+j)
                stTime(iRegion+1+j) = NaN; 
                endTime(iRegion+j) = NaN;
                j = j + 1;
            end
            if iRegion+1+j > length(stTime)
                break
            end
%             IdxList = cellfun(@transpose,CC.PixelIdxList,'UniformOutput',false);
%             certain = GreyImage([IdxList{iRegion:iRegion+j}]);
%             for n = iRegion:iRegion+j
%                 [I,~] = ind2sub(size(GreyImage),CC.PixelIdxList{n}(find(GreyImage(CC.PixelIdxList{n}) == max(certain))));
%                 if ~isempty(I)
%                     break
%                 end
%             end
           instantBinary = imageIn(:,Boundary(iRegion,1):Boundary(iRegion+j,1)+Boundary(iRegion+j,3));
           instantPattern = GreyImage(:,Boundary(iRegion,1):Boundary(iRegion+j,1)+Boundary(iRegion+j,3));
        end
        iRegion = iRegion + j + 1;
        j = 0;
%         cFreq = [cFreq freq(I(1))];
        S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});
        centroids = cat(1, S.WeightedCentroid);
        cFreq = [cFreq freq(round(centroids(2)))];
    end
    stTime = stTime(~isnan(stTime))';
    endTime = endTime(~isnan(endTime))';
    stTime = stTime(2:end-1);
    endTime = endTime(2:end-1);
    %         freqLower = min([length(pseudoFreq)*ones(length(validRegionIdx),1),...
    %             round( centroids(validRegionIdx(:),2)+boundary(validRegionIdx(:),4)/2)]') ;
    %         freqUpper = max([1*ones(length(validRegionIdx),1),...
    %             round( centroids(validRegionIdx(:),2)-boundary(validRegionIdx(:),4)/2)]') ;
    %         bwFreq = pseudoFreq(freqUpper ) - pseudoFreq(freqLower) ;
    Dfre{iChannel} = cFreq;
    Duration{iChannel} =  endTime-stTime ; % second
    Start = stTime;
    End =  endTime;
    interval{iChannel} = Start(2:end) - End(1:end-1);

    i = i + 1;
    toc
end
%%

% close all
figure;
subplot(1,3,1)
histogram(Dfre,50,'normalization','Probability')
xlabel('Drift Frequency(Hz)')
ylabel('Probability')
subplot(1,3,2)
histogram(Duration,50,'normalization','Probability')
xlabel('Duration(ms)')
ylabel('Probability')
subplot(1,3,3)
histogram(interval,50,'normalization','Probability')
xlabel('Interval(ms)')
ylabel('Probability')

%%
figure;
subplot(1,3,1)
[x,n] = histcounts(Dfre,50,'normalization','Probability')
loglog(n(1:end-1),x,'o')
xlabel('Drift Frequency(Hz)')
ylabel('Probability')
subplot(1,3,2)
[x,n] = histcounts(Duration,50,'normalization','Probability')
loglog(n(1:end-1),x,'o')
xlabel('Duration(ms)')
ylabel('Probability')
subplot(1,3,3)
[x,n] = histcounts(interval,50,'normalization','Probability')
loglog(n(1:end-1),x,'o')
xlabel('Interval(ms)')
ylabel('Probability')


% figure;
f = fitdist(Dfre','stable') 
f = fitdist(Duration','stable') 
f = fitdist(interval','stable') 


%% for bandpass
figure;
subplot(2,2,1)
histogram(GammaBurstEvent.burst_du_steps{1},50,'normalization','Probability')
xlabel('Duration(s)')
ylabel('Probability')
subplot(2,2,3)
histogram(GammaBurstEvent.flat_du_steps{1},50,'normalization','Probability')
xlabel('Interval(s)')
ylabel('Probability')

subplot(2,2,2)
[x,n] = histcounts(GammaBurstEvent.burst_du_steps{1},50,'normalization','Probability')
loglog(n(1:end-1),x,'o')
xlabel('Duration(s)')
ylabel('Probability')
subplot(2,2,4)
[x,n] = histcounts(GammaBurstEvent.flat_du_steps{1},50,'normalization','Probability')
loglog(n(1:end-1),x,'o')
xlabel('Interval(s)')
ylabel('Probability')


f = fitdist(GammaBurstEvent.burst_du_steps{1}','stable') 
f = fitdist(GammaBurstEvent.flat_du_steps{1}','stable') 


%% multitaper method
movingwin=[0.2 0.002]; % set the moving window dimensions 30ms*6
params.Fs=fsTemporal; % sampling frequency
params.fpass=[30 80]; % frequencies of interest
params.tapers=[3 3]; % tapers
params.trialave=1; % average over trials
params.err=0; % no error computation
data=sigIn; % data from channel 1
[S1,t,f]=mtspecgramc(data,movingwin,params);

%%
figure;
plot_matrix(S1(1000:2000,:),t(1000:2000),f);
xlabel([]); % plot spectrogram
% caxis([8 28]); 
colorbar;

%%
DurMat = cell2mat(Duration) ;
interMat = cell2mat(interval) ;
DuFull = [DurMat, -DurMat] ;
interFull = [interMat, -interMat] ;
%% fitting
sigDis = (patternScale) ;
sigDis(sigDis<=0) = [];
sigFull = [sigDis,-sigDis] ;
close all
figure;
[x,n] = histcounts(sigDis,20,'normalization','pdf') ;
semilogy((n(1:end-1)+n(2:end))/2,x,'.','markersize',12,'linewidth',2);
f1 = fitdist(sigFull','stable')
x1 = 1:2:100 ;
y1 = 2*pdf(f1,x1) ;
hold on
semilogy(x1,y1,'--','linewidth',1)
f2 = fitdist(sigDis','exponential');
y2 = pdf(f2,x1) ;
hold on
semilogy(x1,y2,'*-','linewidth',1)
f3 = fitdist(sigDis','lognormal');
y3 = pdf(f3,x1) ;
hold on
semilogy(x1,y3,'o-','linewidth',1)

% negloglik(y1)
% wbllike(y1.ParameterValues,x1)
legend('original data','alpha stable','exponential','log-normal')
%legend('original data','alpha stable','exponential')

%%
figure;
% hold on
[x,n] = histcounts(sigDis,200,'normalization','pdf') ;
loglog((n(1:end-1)+n(2:end))/2,x,'.','markersize',12,'linewidth',2);