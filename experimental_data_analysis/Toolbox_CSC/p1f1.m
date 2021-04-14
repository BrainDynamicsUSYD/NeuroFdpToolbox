function p1f1
%% for project 1, figure 1. 
%% load data
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

%% Find Gamma burst in the time domain for fig.1
numChannels = size(sigIn,1)*size(sigIn,2) ;   % sigIn is the abs hilbertSig
sigReshape = reshape(sigIn,numChannels,[]) ;
stdVal = 2.5;
GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,0,badChannels,...
    numChannels,stdVal) ;
%% select a row for visualization
selectCol = 5 ;         % 5 for the paper
selectBur = 16180 ;       % 16180 for the paper
selectColTF = 1 ;       % 1 for the paper
numRow = size(sigIn,2) ;
sumRowBurst = sum(GammaBurstEvent.is_burst((selectCol-1)*numRow+1:selectCol...
    *numRow,:),1) ;
burstTimeIdx =  find(sumRowBurst>2) ;    % guarantee multiple bursts (optional)

if length(burstTimeIdx) < selectBur
    error('there is not that many bursts! Reduce selectBur.')
end
% back to sigOri
burstSelect = burstTimeIdx(selectBur)+fix(fsTemporal) + fix(0.2*fsTemporal) ;   

burstIdx = find(GammaBurstEvent.is_burst((selectCol-1)*numRow+1:...
    selectCol*numRow,burstTimeIdx(selectBur)) == 1) ;
if length(burstIdx) < selectColTF
    error('there is not that many bursts! Reduce selectColTF.')
end
burstTFSelect = burstIdx(selectColTF) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig.1 in the paper is generated in this section
close all
% fig.1(B) Broadband signals
figure;
set(gcf,'Position',[675 541 1006 420])
halfWidth = fix(0.3*fsTemporal) ;
plotStart = burstSelect-halfWidth ;
plotEnd = burstSelect+halfWidth ;

count = 1 ;
for xChan = 1:numRow
    for yChan = selectCol
        plot(plotStart:plotEnd,squeeze(sigOri(xChan,yChan,plotStart:plotEnd)) - ...
            0.5*(count-1)*max(sigOri(:)),'linewidth',1)
        hold on
        legendName{count} = ['(',num2str(xChan),', ',num2str(yChan),')'] ;
        
        count = count + 1 ;
    end
end
% legend(legendName)
% title('time series of the raw signals over 20s (spacing 2 electrode)')
% grid on
% xlabel('time (ms)')
% ylabel('normalised amplitude')
axis off ; 

% fig.1(C) Gamma bandpass signals
figure;
set(gcf,'Position',[675 541 1006 420])
plotStartG = plotStart-fix(fsTemporal)+1 ;   % shift for bandpass signals
plotEndG = plotEnd-fix(fsTemporal)+1 ;

count = 1 ;
for xChan = 1:numRow
    for yChan = selectCol
        plot(plotStart:plotEnd,squeeze(bandpassSig(1,xChan,yChan,plotStartG:plotEndG)) - ...
            0.5*(count-1)*max(bandpassSig(:)),'linewidth',1)
        hold on
        legendName{count} = ['(',num2str(xChan),', ',num2str(yChan),')'] ;
        
        count = count + 1 ;
    end
end
% figure;
% plotStart = burstSelect+fix(fsTemporal)+fix(0.2*fsTemporal)+1-fix(3*fsTemporal) ;
% plotEnd = burstSelect+fix(fsTemporal)+fix(0.2*fsTemporal)+1+fix(3*fsTemporal) ;
% 
% count = 1 ;
% for xChan = 1:10
%     for yChan = 4
%         plot(squeeze(sigIn(xChan,yChan,plotStart:plotEnd)) - (count-1)*max(sigIn(:)))
%         hold on
%         legendName{count} = ['(',num2str(xChan),', ',num2str(yChan),')'] ;
%         
%         count = count + 1 ;
%     end
% end
% legend(legendName)
% title('time series of the raw signals over 20s (spacing 2 electrode)')
% grid on
% xlabel('time (ms)')
% ylabel('normalised amplitude')
axis off ; 

% % fig.1(D) Time-space analysis
% figure;
% numChannelS = size(sigSmooth,1) ;
% resizeScale = size(sigSmooth,1)/size(sigOri,1) ;
% timeStartW = plotStart-fix(fsTemporal)  ;
% timeEndW = plotEnd+fix(fsTemporal)   ;
% set(gcf,'Position',[675 541 1006 420])
% plotStartTS = burstSelect-fix(fsTemporal) - fix(0.2*fsTemporal)+1-halfWidth ;
% plotEndTS = burstSelect-fix(fsTemporal) - fix(0.2*fsTemporal)+1+halfWidth ;
% timePlot = fix(fsTemporal) + 1 : timeEndW - timeStartW - fix(fsTemporal) ; 
% tempTimeAxis = linspace(plotStart,plotEnd,length(timePlot)) ;
% % imagesc(tempTimeAxis,1:numCol,squeeze(sigIn(1:numCol,selectRow,plotStartTS:plotEndTS)))
% imagesc(tempTimeAxis,1:numRow,squeeze(sigSmooth(1:numChannelS,selectCol*...
%     resizeScale,plotStartTS:plotEndTS)))
% ylabel('Space (electrode)')
% % xlabel('Time (s)')
% % set(gca,'YDir','normal')
% % colorbar
% % axis off

% fig.1(E) Time-frequency analysis (wavelet)
figure;
set(gcf,'Position',[675 541 1006 420])
timeStartW = plotStart-fix(fsTemporal)  ;
timeEndW = plotEnd+fix(fsTemporal)   ;

% freqRange = 30:80 ;
% fc = centfrq('cmor1.5-1') ;
% scalerange = fc./(freqRange/fsTemporal) ;
% scales = scalerange(end):0.5:scalerange(1) ;
% pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/fsTemporal) ;
% 
% tempData =  squeeze (sigOri(burstTFSelect,selectCol,timeStartW:timeEndW)) ;
% 
% wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
% timePlot = fix(fsTemporal) + 1 : timeEndW - timeStartW - fix(fsTemporal) ; 
% tempTimeAxis = linspace(plotStart,plotEnd,length(timePlot)) ;
% uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,timePlot)) )))
% xlim([min(tempTimeAxis),max(tempTimeAxis)])
%c
% % xlabel('Time (s)')
% set(gca,'YDir','normal')
% % colorbar
% % axis off
% % tidyPlotForIllustrator
% % figs_to_ps_2015

tempData = squeeze (sigOri(burstTFSelect,selectCol,:)) ;
[wt2,f22,coi] = cwt(tempData,fsTemporal,'VoicesPerOctave',30);
    for j = 1:length(coi)
        ind = find(f22<=coi(j));
        wt2(ind,j) = NaN;
    end
%     wt2(f22<1 | f22>100,:) = [] ;
wt2(f22<30 | f22>80,:) = [] ;
    f2 = f22 ;
%     f2(f22<1 | f22>100,:) = [] ;
f2(f22<30 | f22>80,:) = [] ;
    
% close all
timeStartW = plotStart-fix(fsTemporal)  ;
timeEndW = plotEnd+fix(fsTemporal)   ;
timePlot = timeStartW: timeEndW  ; 
tempTimeAxis = linspace(plotStart,plotEnd,length(timePlot)) ;
figure;
imagesc(tempTimeAxis,flip(f2),((abs(wt2(end:-1:1,timeStartW:timeEndW)) )))
set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colorbar
 

%% for submission to PNAS
close all
figure_width = 8.7; %cm  11.4 or 17.8 for two columns
figure_hight = 15; %cm
figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

subplot(3,1,1)
count = 1 ;
for xChan = 4:7
    for yChan = selectCol
        plot(plotStart:plotEnd,squeeze(sigOri(xChan,yChan,plotStart:plotEnd)) - ...
            0.5*(count-1)*max(sigOri(:)),'linewidth',1)
        hold on
        legendName{count} = ['(',num2str(xChan),', ',num2str(yChan),')'] ;
        
        count = count + 1 ;
    end
end
xlim([plotStart,plotEnd])
% axis off ; 

subplot(3,1,2)
count = 1 ;
for xChan = 4:7
    for yChan = selectCol
        plot(plotStart:plotEnd,squeeze(bandpassSig(1,xChan,yChan,plotStartG:plotEndG)) - ...
            0.5*(count-1)*max(bandpassSig(:)),'linewidth',1)
        hold on
        legendName{count} = ['(',num2str(xChan),', ',num2str(yChan),')'] ;
        
        count = count + 1 ;
    end
end
xlim([plotStart,plotEnd])
% axis off ; 

subplot(3,1,3)
addpath([pwd,'/ToolOthers/uimage'])
uimagesc(tempTimeAxis,flip(f2),((abs(wt2(end:-1:1,timeStartW:timeEndW)) )))
set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colormap(jet)
 % colorbar

% plot(1:10)
% text(2.5,2,{'The figure size should',' be 8.4 by 15 cm'})
set(gca,'FontSize', 8);  % 6 points for x-axis tickmark labels
xlabel('30 point label', 'fontsize', 8 ); % this must be after the above line!

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc figure_size_control % this is the trick!! 
    
