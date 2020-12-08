%% for project 1, figure 2. 
%% load data, with spatial smoothed 8
addpath(genpath([pwd,'/Toolbox_CSC']))         % Tool
addpath(genpath([pwd,'/ToolOthers/nanconv']))

dataFileName = 'ma027_032' ;
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
surSig = generateSur(sigOri,0,badChannels) ;
% apply bandpass filtering for Gamma signals
subBand = [30,80] ;
bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels) ;
% Hilbert transform for analytic signals
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
% find amplitdue of the analytic Gamma as the input
sigIn = abs(squeeze((hilbertSig))) ;

resizeScale = 8 ;
% sigSmooth = sigIn ;

% resizeScale = 2 ;
numChannelS = resizeScale * size(sigIn,1);
sigSmooth = zeros(numChannelS,numChannelS,size(sigIn,3)) ;

for iTime = 1:size(sigIn,3)
    %Try smoothing
    filtWidth = 3;
    filtSigma = 0.6;   % 0.6
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(:,:,iTime))),imageFilter,'edge', 'nonanout');
    sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);
end
clearvars smoothTemp

%%
figure_width = 17; %cm  11.4 or 17.8 for two columns
figure_hight = 13; %cm
figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

% load the burst statistics (8 times resize for visualisation)
load([pwd,'/Results_data/Project1/SmoothedGamma/v4_00.6SmRes8_2.5SD_ma027_032_Band30HzOctober01_23:56.mat'])

%
figure;
set(gcf,'Position',[3 198 1364 471])
subplot(2,4,1)
flagType = 1 ;
movie_burst_2(sigSmooth,WCentroids,patternIdx,rangeFrameTemp,flagType)
subplot(2,4,2)
flagType = 2 ;
movie_burst_2(sigSmooth,WCentroids,patternIdx,rangeFrameTemp,flagType)
subplot(2,4,3)
% load the burst statistics (2 times resize)
load([pwd,'/Results/Project1/project1_results_recent/final_2.5SD/v4_00.6Sm_2.5SD_ma027_032_Band30HzSeptember23_11:07.mat'])

sigDist = Duration ;
numPts = 14 ;
stableFit(sigDist,numPts)
subplot(2,4,4)
sigDist = centInterval ;
numPts = 14 ;
stableFit(sigDist,numPts)

subplot(2,3,4)
p1f2_msdSample(WCentroids)
subplot(2,3,5)
[superIdx,meanSlope] = p1f2_msdDist(WCentroids) ;
subplot(2,3,6)
p1f2_msdCollapse(WCentroids,superIdx)

set(gca,'FontSize',10)
set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc figure_size_control % this is the trick!! 
