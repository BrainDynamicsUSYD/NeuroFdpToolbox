function job_different_freq_threshold(arrayID)

cd ..
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%
% arrayID = 8 ;

switch arrayID
    case 1
        dataFileName = 'ma027_032' ;
        % dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 1;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
    case 2
        dataFileName = 'ma027_032' ;
        % dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 2;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
    
    case 3
        dataFileName = 'ma027_032' ;
        % dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 3;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
    case 4
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 4;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
    case 5
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 5;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
    case 6
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 6;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
               
     case 7
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 7;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;        
        
     case 8
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 8;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
     case 9
        dataFileName = 'ma027_032' ;
        dataFileName = 'ma025_03' ;
        % timeStart = 0 ;     % for finding centres
        timeStart2 = 9;
        subBand = [30 , 80] ;
        plotAmp = 1 ;
        % timeStep = 5 ;
        makeMovie = 1;
        makeMovieSmooth = 0;
        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')

%%
if plotAmp
    dataFileName = ['Amp',dataFileName] ;
else
    dataFileName = ['Phase',dataFileName] ;
end

%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
% clearvars LFPs
% sigOriTemp = reshape(sigOri,10*10,[]) ;
% surSigTemp = nan(size(sigOriTemp)) ;

%%
[sigOriRaw,~,~] = preprocess_LFP(LFPs, 0) ;

%% fft
subBand = [0.5,4 ; 4,9 ; 9, 13; 13, 30; 30, 80; 80, 200] ;
subBand = [1,2 ; 3,4 ; 30, 40; 50, 60; 70, 80; 90, 100] ;   % Delta and Gamma
subBand = [1,3 ; 3,5 ; 7, 9; 15, 17; 31, 33; 63, 65] ;      % log spaced
subBand = [1,3 ; 9,11 ; 17, 19; 25, 27; 33, 35; 41, 43] ;   % linear spaced
% subBand = [1,10 ; 11,20 ; 21, 30; 31, 40; 41, 50; 51, 60] ;      % linear spaced/large


nfft = size(sigOri,3) ;
vpf = zeros(nfft,1) ;
for xPoint = 1:10
    for yPoint = 1:10
if(isnan(sigOri(xPoint,yPoint,1)))
    continue
end
% You can compute the Fourier transform of an image in Matlab using fft2 as follows:
vf = fft(squeeze(sigOri(xPoint,yPoint,:)),nfft) /sqrt(nfft);
%vf2 = fft(squeeze(filteredLFPs(xPoint,yPoint,:)),nfft) ;
% The function fftshift is needed to put the DC component (frequency = 0) in the center.  
%The power spectrum is just the square of the modulus of the Fourier transform,  which is obtained as follows:
vpf = vpf+abs(vf).^2 ;

    end
end
vpf = vpf/(100-length(badChannels)) ;
%vpf2 = abs(vf2).^2 ;
% To display the power spectrum you will need to define some frequency coordinates.  
% If the image is of size N, then the frequencies run from -N/2 to N/2-1 (assuming N is even):
% [~,~,N] = size(LFPs) ;
freqIndex = 1 : ceil(nfft/2) ;
freqAxis = freqIndex/nfft*Fs ;
% Then display the log power spectrum using imagesc:

figure; set(gcf,'Visible','On')
% plot(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
loglog(freqAxis,vpf(freqIndex),'o'); 
xlim([0.1 400])
xlabel('frequency (Hz)')
ylabel('power')
title('power spectrum')

hold on
%figure; set(gcf,'Visible','On')
% plot(freqAxis,vpf(freqIndex))
%subplot(1,2,1)
%loglog(freqAxis,vpf2(freqIndex)); 
%xlim([1 500])

p = polyfit(log(freqAxis),log(vpf(freqIndex)'),1) ;
    y = exp(polyval(p,log(freqAxis))) ;
    slope = p(1) ;
    loglog(freqAxis,y)

plotMax = zeros(size(subBand,1)-1,1) ;
freqAverage = [4,8,12,16,30] ;
for iLoop = 1:length(plotMax)
    [~,idxFreq] = min(abs(freqAxis - freqAverage(iLoop))) ;
    plotMax(iLoop) = y(idxFreq) ;
end

%% smooth the original signal
% addpath(genpath([pwd,'/ToolNeuroPatt']))
% addpath(genpath([pwd,'/ToolOthers/nanconv']))
% sigOriSmooth = [] ;
% resizeScale = 2 ;
% for iTime = 1:fix(80*fsTemporal) %size(sigOri,3)% plotLength
%     timeSlot = iTime+timeStart ;
%      %Try smoothing
%     filtWidth = 3;
%     filtSigma = 0.6;
%     imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%     smoothTemp = nanconv(abs(squeeze(sigOri(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
%     sigOriSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);   
% end

%%

% totalSur = 1 ;
% for nSur = 1:totalSur
%     randNumHalf =  2*pi*rand(size(sigOri,3)/2-1,1) ;
%     randNum = [0;randNumHalf;0;-flip(randNumHalf) ]' ;
%     for iChannel = setdiff(1:100,badChannels)
%         freqSig = fft(sigOriTemp(iChannel,:)) ;
%         absFreq = abs(freqSig) ;
%         phaseFreq = angle(freqSig) ;
%         reconSig = ifft(absFreq.*exp(1i*phaseFreq)) ;
%         
%         surSigTemp(iChannel,:) = ifft(absFreq.*exp(1i*(phaseFreq+randNum))) ;
%     end
%  surSig = reshape(real(surSigTemp),10,10,[]) ;

% surSig = sigOri ;
% find subband signals
bandpassSig = find_bandpassSig(sigOri,subBand, fsTemporal,3,badChannels,0) ;
%%
discardTimeSteps = fix(1*fsTemporal) ;
outSection = discardTimeSteps+1 : (size(sigOri,3)-discardTimeSteps);
sigOriMatchSize = reshape(sigOri(:,:,outSection),1,10,10,[]) ;
%
bandpassSig = [sigOriMatchSize;bandpassSig] ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4,badChannels) ;

%%
if plotAmp
    sigIn = abs(squeeze(hilbertSig(:,:,:,fix(11*fsTemporal)+1:fix(50*fsTemporal)) ));
else
    sigIn = angle(squeeze(hilbertSig(:,:,:,fix(11*fsTemporal)+1:fix(50*fsTemporal)) ));
end
sigOriStart = fix(1 * fsTemporal) + fix(0.2 * fsTemporal) + fix(11*fsTemporal) ;

% clearvars bandpassSig hilbertSig


%% Interpolation and smoothing

close all
% plotLength = fix(100*fsTemporal) ;
sigSmooth = zeros(size(sigIn,1),20,20,size(sigIn,4)) ;
% timeStart = 5000 ;

addpath(genpath([pwd,'/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))

resizeScale = 2 ;
for bandIdx = 1:size(sigIn,1)
for iTime = 1:size(sigIn,4)% plotLength
    timeSlot = iTime ;
     %Try smoothing
    filtWidth = 3;
    filtSigma = 0.6;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(bandIdx,:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    sigSmooth(bandIdx,:,:,iTime) = imresize(smoothTemp, resizeScale);   
end
end
clearvars smoothTemp


%% find 95 percentile
sigPlot = zeros(size(sigSmooth)) ;
for bandIdx = 1:size(sigSmooth,1)
    sigSmoothTemp = sigSmooth(bandIdx,:,:,:) ;

if ~plotAmp
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>pi*5/6) = 1 ;
    sigPlot(bandIdx,:,:,:) = sigBinary.*sigSmoothTemp ;
else
    boundPrctile = 95 ;
    prcBound(bandIdx,:) = prctile(sigSmoothTemp(:),boundPrctile) ;
    sigBinary = zeros(size(sigSmoothTemp)) ;
    sigBinary(sigSmoothTemp>prcBound(bandIdx,:)) = 1;
    sigPlot(bandIdx,:,:,:) = sigBinary.*sigSmoothTemp ;
end
end
%%
%     CC = bwconncomp(sigBinary) ;                % 3D
%     B = regionprops(CC,'BoundingBox');
%     boundary = cat(1, B.BoundingBox);
%     Area = regionprops(CC,'Area') ;
%     count = 1;                                 % for counting bursts
%     Centroids = [] ;
%     WCentroids = [] ;
%     instantScale = [] ;
%     instantPeakAmp = [] ;
%     instantTotalPower = [] ;
%     rangeFrame = [] ;
%     % Amp = [] ;
%     
%     burstIdxOld = zeros(size(sigPlot)) ;
%     WCentroidsPlot = cell(size(sigBinary,3),1) ;
%     %burstPlot = zeros(size(sigIn)) ;
%     
%     for iBurst = 1: size(CC.PixelIdxList,2)
%         currentIdx = CC.PixelIdxList{iBurst} ;
%         % Duration1(iBurst) = boundary(iBurst,6) ;    % should be similar to
%         % Duration2
%         burstTimeEnd = floor((currentIdx(end)-1)/100) +1 ;
%         burstTimeStart = floor((currentIdx(1)-1)/100) +1 ;
%         Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;    %
%         if Duration2(iBurst) < minBurstTime
%             continue
%         end
%         % Amp = [Amp; sigPlot(currentIdx)];
%         
%         Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration
%         
%         patternScale(count) = Area(iBurst).Area ;  % 4 total scale
%         
%         burstIdxTemp = zeros(size(sigPlot)) ;
%         burstIdxTemp(currentIdx) = 1 ;
%         currentBurst = sigPlot.*burstIdxTemp ;
%         sumAmp(count) = sum(currentBurst(:)) ;     % sum of amplitude
%         peakAmp(count) = max(currentBurst(:)) ;    % peak amplitude
%       % 
%         burstIdx = burstIdxTemp+burstIdxOld ;
%         burstIdxOld = burstIdx;
%         
%         timeCount = 1 ;
%         
%         for iTime = burstTimeStart:burstTimeEnd  
%             instantPattern = currentBurst(:,:,iTime) ; 
%             % burstPlot(:,:,iTime) = burstPlot(:,:,iTime)+instantPattern ;
%             instantBinary = burstIdxTemp(:,:,iTime) ;
%             instantScale{count}(timeCount,:) = sum(instantBinary(:)) ;  % instant scale
%             instantPeakAmp{count}(timeCount,:) = max(max(currentBurst(:,:,iTime))) ;
%             instantTotalPower{count}(timeCount,:) = sum(sum(currentBurst(:,:,iTime).^2) ) ;
%             S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});
%             Centroids{count}(timeCount,:) = cat(1, S.Centroid);
%             WCentroids{count}(timeCount,:) = cat(1, S.WeightedCentroid) ;
%             timeCount = timeCount + 1;
%             
%             WCentroidsPlot{iTime} = [WCentroidsPlot{iTime},cat(1, S.WeightedCentroid)] ;
%         end
%         rangeFrame(count,:) = [burstTimeStart,burstTimeEnd] ; 
%         if count>1
%             firstCentroidsLoc = squeeze(WCentroids{count}(1,:)) ;
%             distCent(count) = sqrt(sum((firstCentroidsLoc - lastCentroidsLoc).^2)) ;
%             
%             firstCentroidsTime = burstTimeStart ;
%             centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
%         end
%         lastCentroidsLoc = squeeze(WCentroids{count}(end,:)) ;
%         lastCentroidsTime = burstTimeEnd ;
%         count = count+1 ;
%         
%     end


%% make a movie
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;

% sigBurst = sigIn.*burstIdx ;
% timeStart2 = fix(fsTemporal) ;
if makeMovie
    vidTitle = [pwd,'/Results/Project2/Results/moviesRaw/',dataFileName,'_',...
        num2str(timeStart2),'0s_',dateStr,'.avi'] ;
    vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 20 ;
    open(vidObj);
     fig=figure ;
     set(gcf,'Position',[255 35 1791 919])     
    
    for iTime = 1:10000
        timeSlot = iTime+fix(timeStart2*fsTemporal*10) ;
        % delta
        subplot(2,3,1)
        imagesc(squeeze(sigPlot(2, :,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Delta band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        title(['1-2Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear
        title(['1-3Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % log
        % title(['1-3Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear,large
        % title(['1-10Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        
        if plotAmp
            % tempSig = sigPlot(2, :,:,:) ;
            % caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(2,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,3,2)
        imagesc(squeeze(sigPlot(3, :,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Theta band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        title(['3-4Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear
        title(['9-11Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % log
        % title(['3-5Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear,large
        % title(['11-20Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        
        if plotAmp
%             tempSig = sigPlot(3, :,:,:) ;
%             caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(3,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,3,3)
        imagesc(squeeze(sigPlot(4, :,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Alpha band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        title(['30-40Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear
        title(['17-19Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % log
        %title(['7-9Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear,large
        % title(['21-30Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        
        if plotAmp
%             tempSig = sigPlot(4, :,:,:) ;
%             caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(4,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,3,4)
        imagesc(squeeze(sigPlot(5, :,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Beta band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        title(['50-60Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear
        title(['25-27Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % log
        % title(['15-17Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear,large
        % title(['31-40Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        
        if plotAmp
%             tempSig = sigPlot(5, :,:,:) ;
%             caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(5,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        
        subplot(2,3,5)
        imagesc(squeeze(sigPlot(6,:,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Gamma band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        title(['70-80Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear
        title(['33-35Hz band at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        % log
        %title(['31-33Hz band at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
        % linear,large
        % title(['41-50Hz band at ', ...
         %   int2str(timeSlot/fsTemporal*1000),'ms'])
        
        if plotAmp
%             tempSig = sigPlot(6, :,:,:) ;
%             caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(6,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,3,6)
        % imagesc(sigOri( :,:,sigOriStart+iTime+timeStart))
        imagesc(squeeze(sigPlot(1,:,:,timeSlot)))
        set(gca,'YDir','normal')
        title(['Raw signal at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        if plotAmp
%             tempSig = sigPlot(1, :,:,:) ;
%             caxis([0 0.5*max(tempSig(:))])
            caxis([0 2*prcBound(1,:)])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        % pause(0.01)
        cla
        
    end
    close(vidObj);
end

%% make a movie
% close all
% % sigBurst = sigIn.*burstIdx ;
% if makeMovieSmooth
%     t = datetime('now') ;
%     dateStr = datestr(t,'mmmmdd_HH:MM') ;
% 
%     dataFileName = ['Smooth',dataFileName] ;
%     vidTitle = [pwd,'/Results/movies/New',dataFileName,'_',...
%         num2str(timeStart2),dateStr,'.avi'] ;
%     vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);
%      fig=figure ;
%      set(gcf,'Position',[255 35 1791 919])     
%     
%     for iTime = 1:1000
%         timeSlot = iTime+timeStart2 ;
%         % delta
%         subplot(2,3,1)
%         imagesc(squeeze(abs( sigSmooth(1, :,:,iTime+timeStart) ).^2))
%         set(gca,'YDir','normal')
%         title(['Delta band at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             tempSig = sigSmooth(1, :,:,:) ;
%             caxis([0 plotMax(1)/8])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         subplot(2,3,2)
%         imagesc(squeeze(abs(sigSmooth(2, :,:,iTime+timeStart) ).^2))
%         set(gca,'YDir','normal')
%         title(['Theta band at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             tempSig = sigSmooth(2, :,:,:) ;
%             caxis([0 plotMax(2)/8])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         subplot(2,3,3)
%         imagesc(squeeze(abs(sigSmooth(3, :,:,iTime+timeStart) ).^2))
%         set(gca,'YDir','normal')
%         title(['Alpha band at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             tempSig = sigSmooth(3, :,:,:) ;
%             caxis([0 plotMax(3)/8])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         subplot(2,3,4)
%         imagesc(squeeze(abs(sigSmooth(4, :,:,iTime+timeStart)) .^2))
%         set(gca,'YDir','normal')
%         title(['Beta band at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             tempSig = sigSmooth(4, :,:,:) ;
%             caxis([0 plotMax(4)/8])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         
%         subplot(2,3,5)
%         imagesc(squeeze(abs(sigSmooth(5,:,:,iTime+timeStart) ).^2))
%         set(gca,'YDir','normal')
%         title(['Gamma band at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             tempSig = sigSmooth(5, :,:,:) ;
%             caxis([0 plotMax(5)/8])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         subplot(2,3,6)
%         imagesc(abs(sigSmooth( :,:,sigOriStart+iTime+timeStart)).^2)
%         set(gca,'YDir','normal')
%         title(['Raw signal at ', ...
%             int2str(timeSlot/fsTemporal*1000),'ms'])
%         if plotAmp
%             caxis([0 max(sigSmooth(:).^2/8 )])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
%         
%         
%         writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
%         pause(0.01)
%         cla
%         
%     end
%     close(vidObj);
% end