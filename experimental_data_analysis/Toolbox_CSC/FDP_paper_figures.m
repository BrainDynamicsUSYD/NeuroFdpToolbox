% This is the main function to generate all the figures in the paper:
% Lévy walk dynamics of gamma burst patterns in primate cerebral cortex
% submitted to Communication Biology
%
% Author: Xian Long; xlon3884@uni.sydney.edu.au
%
addpath(genpath([pwd,'/Toolbox_CSC']))         % add sub-functions
addpath(genpath([pwd,'/ToolOthers/nanconv']))  % to interpolate missing chans

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 1. schematic diagram of Levy flight
schematicLF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 2. gamma bursts
p1f1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig 3A. pattern trajectories
% flagDetectPatt = 1 ;         % 1 run the pattern detection, 0 load directly
% if flagDetectPatt
%     dataFileName = 'ma027_032' ;
%     stdPara = 2.5 ;
%     [Duration,patternScale,centInterval,WCentroids] = 
% else
%     
% end

%% fig 3B. pattern statistics
DurationFull = [] ;
flagDetectPatt = 0 ;         % 1 run the pattern detection, 0 load directly
if flagDetectPatt
    Project1(1)
    Project1(2)
    Project1(3)
    Project1(4)
else
    load([pwd,'/Results_data/Project1/Project1_submit/',...
        '00.6SmRes2_2.5SD_my147_53_Band30HzMay11_22:53.mat_v4.mat'])
    DurationFull = [DurationFull,Duration] ;
    load([pwd,'/Results_data/Project1/Project1_submit/',...
        '00.6SmRes2_2.5SD_ma025_03_Band30HzMay11_22:49.mat_v4.mat'])
    DurationFull = [DurationFull,Duration] ;
    load([pwd,'/Results_data/Project1/Project1_submit/',...
        '00.6SmRes2_2.5SD_my144_101_Band30HzMay11_12:01.mat_v4.mat'])
    DurationFull = [DurationFull,Duration] ;
    load([pwd,'/Results_data/Project1/Project1_submit/',...
        '00.6SmRes2_2.5SD_ma027_032_Band30HzMay11_22:37.mat_v4.mat'])
    DurationFull = [DurationFull,Duration] ;
    sigDist = DurationFull ;

    figure;
    subplot(2,2,1)
    stableFit(sigDist,numPts,1)
    subplot(2,2,2)
    stableFitOne(sigDist,numPts,1)
    subplot(2,2,3)
    stableFitCDFTwo(sigDist,numPts,1)
    subplot(2,2,4)
    stableFitCDFOne(sigDist,numPts,1)
    
    pf_power = @(x,alpha,xmin) x.^-alpha*(alpha-1)*xmin^(alpha-1);
    cdf_power = @(x,alpha,xmin) (x/xmin).^(-alpha+1) ;
    [lambdaHat1,lambdaCI] = mle(sigDist, 'pdf',pf_power, 'start',[1,30], 'lowerbound',[1,30], 'upperbound', [5,200]) 
    [x,n] = histcounts(sigDist,100,'normalization','cdf') ;
    n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ;
    sigDist(sigDist<58) = [] ;
    xmin = min(sigDist) ;
    xmax = max(sigDist) ;
    pf_TP = @(x,alpha) (alpha+1)*(xmax.^(alpha+1)-xmin.^(alpha+1)).^(-1).*x.^(alpha) ;
    cdf_TP = @(x,alpha) (x.^(alpha+1)-xmin.^(alpha+1))./(xmax.^(alpha+1)-xmin.^(alpha+1)) ;
    [lambdaHat1,lambdaCI] = mle(sigDist, 'pdf',pf_TP, 'start',[-4], 'lowerbound',[-5], 'upperbound', [-1]) 

    
%  pf_powerC = @(x,alpha,lambda,xmin) x.^-alpha.*exp(-lambda*x)*...
%      lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
% [lambdaHat2,lambdaCI] = mle(sigDist, 'pdf',pf_powerC, 'start',[1,0,30], 'lowerbound',[1,0,30], 'upperbound', [5,10,1000]) ;

%     [x,n] = histcounts(sigDist,100,'normalization','pdf') ;
%     n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ; 
%     figure;
%     loglog(n0,x,'b.','MarkerSize',12)
%     hold on
%     loglog(n0,pf_powerC(n0,lambdaHat2(1),lambdaHat2(2),lambdaHat2(3)),'lineWidth',2)
%     
%     figure;
%     [x,n] = histcounts(sigDist,500,'normalization','pdf') ;
%     n0 = 0.5*((n(2:end))+(n(1:end-1))) ;
%     loglog(n0,x,'b.','MarkerSize',12)
%     hold on
%     loglog(n0,pf_TP(n0,lambdaHat1(1)),'lineWidth',2)
    
    figure;
    [x,n] = histcounts(sigDist,1000,'normalization','cdf') ;
    % n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ;
    n0 = 0.5*((n(2:end))+(n(1:end-1))) ;
    loglog(n0,1-x,'b.','MarkerSize',12)
    hold on
    loglog(n0,1-cdf_TP(n0,lambdaHat1(1)),'r','lineWidth',2)

    [lambdaHat1,lambdaCI] = mle(sigDist, 'pdf',pf_TP, 'start',[-4], 'lowerbound',[-5], 'upperbound', [-1]) 

    paraDist2 = fitdist(sigDist','normal')  ;
    y02 = cdf(paraDist2,n0) ;
    hold on
    loglog(n0,1-y02,'k','lineWidth',2)

        
%     load([pwd,'/Results/Project1/project1_results_recent/final_2.5SD/',...
%         'v4_00.6Sm_2.5SD_ma027_032_Band30HzSeptember23_11:07.mat'])
%         figure;
%     subplot(2,2,1)
%     sigDist = Duration ;
%     stableFit(sigDist,numPts,1)
%     subplot(2,2,2)
%     stableFitOne(sigDist,numPts,1)
%     subplot(2,2,3)
%     stableFitCDFTwo(sigDist,numPts,1)
%     subplot(2,2,4)
%     stableFitCDFOne(sigDist,numPts,1)
end

%% fig 3CDE. MSD
close all
superFull = [];

sdPara = 2;
tic

         load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_my147_53_Band30HzMay11_22:53.mat_v4.mat'],'WCentroids','centInterval')

% step length fitting 1D
% figure;
superFull = [] ;
minDisp = 0.1 ;
xmax = 12.8*200 ; 
noiseThre = 20 ;
grid2Micro = 200 ;
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])
set(gcf,'Position',[534 721 293 245])

%%
         load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_ma025_03_Band30HzMay11_22:49.mat_v4.mat'],'WCentroids','centInterval')

% step length fitting 1D
% figure;
superFull = [] ;
xmax = 12.4*200 ; 
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])
set(gcf,'Position',[534 721 293 245])


%%
         load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_my144_101_Band30HzMay11_12:01.mat_v4.mat'],'WCentroids','centInterval')
                 WCenFull = [WCenFull,WCentroids] ;
        InterFull = [InterFull,centInterval] ;
% step length fitting 1D
% figure;
superFull = [] ;
xmax = 12.0*200 ; 
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])
set(gcf,'Position',[534 721 293 245])

%%
         load([pwd,'/Results_data/Project1/Project1_submit/',...
      '00.6SmRes2_2.5SD_ma027_032_Band30HzMay11_22:37.mat_v4.mat'],'WCentroids','centInterval')
        WCenFull = [WCenFull,WCentroids] ;
        InterFull = [InterFull,centInterval] ;

% step length fitting 1D
% figure;
superFull = [] ;
xmax = 10.0*200 ; 
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])
set(gcf,'Position',[534 721 293 245])



% [superIdx,meanSlope,stdSlope] = p1f2_msdDist(WCentroids) ;
% superFull = [superFull,superIdx] ;



toc
clearvars WCentroids


%% angle model
         load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_ma025_03_Band30HzMay11_22:49.mat_v4.mat'],'WCentroids','centInterval')

% step length fitting 1D
% figure;
superFull = [] ;

iFig = 1 ;
turnAngle  = 20 ;
subplot(2,3,iFig)
% figure;
superIdx = [] ;
flagUseBetween = 1 ;
maxOffset = 200*8 ;   %10  8 for 120
minDisp = 200*0.12 ;
flagPlotPart = 1;

[DisplaceNorm,lambdaHat1] = p1f2_stepLength_TP_2D(WCentroids,superIdx,turnAngle,...
    flagUseBetween,maxOffset,minDisp,flagPlotPart) ;
title(['\theta = ',num2str(turnAngle), char(176) ,' , \lambda = ',num2str(-lambdaHat1,'%0.2f')])

ylim([10^-3 10^0])
xlim([10^1 10^4])

% set(gcf,'Position',[534 721 293 245])
% exportgraphics(gcf,['S2_angle',num2str(turnAngle),'.pdf'],'ContentType','vector')


for turnAngle  = 40:20:120 ;
    iFig = iFig+1 ;
subplot(2,3,iFig)
% figure;
superIdx = [] ;
flagUseBetween = 1 ;
maxOffset = 200*8 ;   %10  8 for 120
minDisp = 200*0.1 ;
flagPlotPart = 1;

[DisplaceNorm,lambdaHat1] = p1f2_stepLength_TP_2D(WCentroids,superIdx,turnAngle,...
    flagUseBetween,maxOffset,minDisp,flagPlotPart) ;
title(['\theta = ',num2str(turnAngle), char(176) ,' , \lambda = ',num2str(-lambdaHat1,'%0.2f')])

ylim([10^-3 10^0])
xlim([10^1 10^4])
% set(gcf,'Position',[534 721 293 245])
% exportgraphics(gcf,['S2_angle',num2str(turnAngle),'.pdf'],'ContentType','vector')


end
annotation(gcf,'textbox',...
    [0.0731343873517787 0.945709281961471 0.0306205533596838 0.0437828371278459],...
    'String','A',...
    'LineStyle','none',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.359695652173913 0.945709281961471 0.0306205533596838 0.0437828371278459],...
    'String','B',...
    'LineStyle','none',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.637363636363635 0.945709281961471 0.0306205533596838 0.0437828371278459],...
    'String','C',...
    'LineStyle','none',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.0731343873517787 0.455341506129597 0.0306205533596838 0.0437828371278459],...
    'String','D',...
    'LineStyle','none',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.359695652173913 0.455341506129597 0.0306205533596838 0.0437828371278459],...
    'String','E',...
    'LineStyle','none',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.637363636363635 0.455341506129597 0.0306205533596838 0.0437828371278459],...
    'String','F',...
    'LineStyle','none',...
    'FitBoxToText','off');

%% s.d.
close all
      load([pwd,'/Project1_preprocessedData/',...
    '00.6SmRes2_1.5SD_ma025_03_Band30HzJanuary16_01:52.mat_v4.mat'],'WCentroids','centInterval')

minDisp = 0.2 ;
xmax = 7.2*200 ; 
noiseThre = 20 ;
grid2Micro = 200 ;
subplot(1,3,1)
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;
ylim([10^-3 10^0])
% set(gcf,'Position',[534 721 293 245])

  load([pwd,'/Project1_preprocessedData/',...
    '00.6SmRes2_2SD_ma025_03_Band30HzJanuary16_01:04.mat_v4.mat'],'WCentroids','centInterval')
minDisp = 0.2 ;
xmax = 10.8*200 ; 
noiseThre = 20 ;
grid2Micro = 200 ;
subplot(1,3,2)
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;
ylim([10^-3 10^0])
% set(gcf,'Position',[534 721 293 245])

     load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_ma025_03_Band30HzMay11_22:49.mat_v4.mat'],'WCentroids','centInterval')

 minDisp = 0.2 ;
xmax = 12.8*200 ; 
noiseThre = 20 ;
grid2Micro = 200 ;
subplot(1,3,3)
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;
ylim([10^-3 10^0])
% set(gcf,'Position',[534 721 293 245])
 
 
%% different threshold
load([pwd,'/Results_data/Project1/Project1_submit/',...
     '00.6SmRes2_2.5SD_ma025_03_Band30HzMay11_22:49.mat_v4.mat'],'WCentroids','centInterval')

      minDisp = 80 ;
xmax = 12*200 ; 
noiseThre = 80 ;
grid2Micro = 200 ;
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])

figure;
      minDisp = 400 ;
xmax = 17*200 ; 
noiseThre = 400 ;
grid2Micro = 200 ;
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])


%% load the data with no spatial scaling
WCenFull = [] ;
sdPara = 1;
tic

figure;
        load([pwd,'/Project1_preprocessedData/',...
    '00.6SmRes1_2.5SD_ma025_03_Band30HzApril06_15:59.mat_v4.mat'],'WCentroids','centInterval')
        WCenFull = [WCenFull,WCentroids] ;

      minDisp = 80 ;
xmax = 6*400 ; 
noiseThre = 80 ;
grid2Micro = 400 ;
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])

figure;
      minDisp = 400 ;
xmax = 17*400 ; 
noiseThre = 400 ;
grid2Micro = 400 ;
DisplaceNorm = p1f2_stepLength_TP(WCentroids,superFull,0,minDisp,xmax,noiseThre,grid2Micro) ;

ylim([10^-3 10^0])
