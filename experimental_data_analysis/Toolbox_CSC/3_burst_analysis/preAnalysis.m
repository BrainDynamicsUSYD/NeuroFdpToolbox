function preAnalysis(sigOri,bandpassSig,sigIn,sigSmooth,badChannels,fsTemporal,stdVal)

%% Find Gamma burst in the time domain for fig.1
numChannels = size(sigIn,1)*size(sigIn,2) ;   % sigIn is the abs hilbertSig
sigReshape = reshape(sigIn,numChannels,[]) ;
% stdVal = 2.5;
GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,0,badChannels,...
    numChannels,stdVal) ;
%% select a row for visualization
selectCol = 5 ;         % 5 for the paper
selectBur = 800 ;       % 800 for the paper
selectColTF = 1 ;       % 5 for the paper
numRow = size(sigIn,2) ;
sumRowBurst = sum(GammaBurstEvent.is_burst((selectCol-1)*numRow+1:selectCol...
    *numRow,:),1) ;
burstTimeIdx =  find(sumRowBurst>2) ;    % guarantee multiple bursts (optional)

if length(burstTimeIdx) < selectBur
    error('there is not that many bursts! Reduce selectBur.')
end
burstSelect = burstTimeIdx(selectBur)+fix(fsTemporal) + fix(0.2*fsTemporal) ;   

burstIdx = find(GammaBurstEvent.is_burst((selectCol-1)*numRow+1:...
    selectCol*numRow,burstSelect) == 1) ;
if length(burstIdx) < selectColTF
    error('there is not that many bursts! Reduce selectColTF.')
end
burstTFSelect = burstIdx(selectColTF) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig.1 in the paper is generated in this section
close all
% fig.1(B) Broadband signals
figure;
set(gcf,'Position',[675 541 1006 420])
halfWidth = fix(0.6*fsTemporal) ;
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
plotStartG = burstSelect-fix(fsTemporal)+1-halfWidth ;
plotEndG = burstSelect-fix(fsTemporal)+1+halfWidth ;

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

% fig.1(D) Time-space analysis
figure;
numChannelS = size(sigSmooth,1) ;
resizeScale = size(sigSmooth,1)/size(sigOri,1) ;
timeStartW = plotStart-fix(fsTemporal)  ;
timeEndW = plotEnd+fix(fsTemporal)   ;
set(gcf,'Position',[675 541 1006 420])
plotStartTS = burstSelect-fix(fsTemporal) - fix(0.2*fsTemporal)+1-halfWidth ;
plotEndTS = burstSelect-fix(fsTemporal) - fix(0.2*fsTemporal)+1+halfWidth ;
timePlot = fix(fsTemporal) + 1 : timeEndW - timeStartW - fix(fsTemporal) ; 
tempTimeAxis = linspace(plotStart,plotEnd,length(timePlot)) ;
% imagesc(tempTimeAxis,1:numCol,squeeze(sigIn(1:numCol,selectRow,plotStartTS:plotEndTS)))
imagesc(tempTimeAxis,1:numRow,squeeze(sigSmooth(1:numChannelS,selectCol*...
    resizeScale,plotStartTS:plotEndTS)))
ylabel('Space (electrode)')
% xlabel('Time (s)')
% set(gca,'YDir','normal')
% colorbar
% axis off

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
    wt2(f22<30 | f22>80,:) = [] ;
    f2 = f22 ;
    f2(f22<30 | f22>80,:) = [] ;
    
% close all
timeStartW = plotStart-fix(fsTemporal)  ;
timeEndW = plotEnd+fix(fsTemporal)   ;
timePlot = timeStartW: timeEndW  ; 
tempTimeAxis = linspace(plotStart,plotEnd,length(timePlot)) ;
figure;
uimagesc(tempTimeAxis,flip(f2),((abs(wt2(end:-1:1,timeStartW:timeEndW)) )))
set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colorbar
 

%% for submission to PNAS
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
% axis off ; 

subplot(3,1,3)
uimagesc(tempTimeAxis,flip(f2),((abs(wt2(end:-1:1,timeStartW:timeEndW)) )))
set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colormap(jet)
 colorbar

% plot(1:10)
% text(2.5,2,{'The figure size should',' be 8.4 by 15 cm'})
set(gca,'FontSize', 8);  % 6 points for x-axis tickmark labels
xlabel('30 point label', 'fontsize', 8 ); % this must be after the above line!

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc figure_size_control % this is the trick!! 
    

%% 
 figure;
     CData2 = abs(wt2);
% CData = abs(wt2);
    Y = prctile(CData2(:),95); % 95
    CData2(CData2 < Y) = 0;
    GreyImage2 = CData2;
    CData2(CData2>0) = 1;
    binaryImage2 = CData2;
    uimagesc(tempTimeAxis,flip(f2),((abs(binaryImage2(end:-1:1,timeStartW:timeEndW)) )))
  set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colorbar  
 
  figure;
     CData2 = abs(wt2);
% CData = abs(wt2);
    % Y = prctile(CData2(:),95); % 95
    Y2 = nanmean(CData2(:)+2*nanstd(CData2(:))) ;
    CData2(CData2 < Y2) = 0;
    GreyImage2 = CData2;
    CData2(CData2>0) = 1;
    binaryImage2 = CData2;
    uimagesc(tempTimeAxis,flip(f2),((abs(binaryImage2(end:-1:1,timeStartW:timeEndW)) )))
  set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colorbar  
 
 
  figure;
     CData2 = abs(wt2);
% CData = abs(wt2);
    % Y = prctile(CData2(:),95); % 95
    Y3 = nanmedian(CData2(:))*6 ;
    CData2(CData2 < Y3) = 0;
    GreyImage2 = CData2;
    CData2(CData2>0) = 1;
    binaryImage2 = CData2;
    uimagesc(tempTimeAxis,flip(f2),((abs(binaryImage2(end:-1:1,timePlot)) )))
  set(gca,'YDir','normal')
 ylabel('Frequency (Hz)')
 xlabel('Time (s)')
 colorbar  
 
%% supplementary figures
% power spectrum around peak
figure;
set(gcf,'Position',[675 541 1006 420])
halfWidth = fix(0.2*fsTemporal) ;
plotStart = burstSelect-halfWidth ;
plotEnd = burstSelect+halfWidth ;
freqAxis = linspace(0,fsTemporal/2,halfWidth+1) ;

powerSp = abs(fft(squeeze(sigOri(xChan,yChan,plotStart:plotEnd)))) ;
loglog(freqAxis,powerSp(1:halfWidth+1))

end