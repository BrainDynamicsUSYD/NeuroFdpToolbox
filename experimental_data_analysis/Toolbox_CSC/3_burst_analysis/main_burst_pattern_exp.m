%% Main Function for Gamma burst pattern
%
% Required sub-folders:
% Toolbox_CSC
% ToolNeuroPatt
% Data
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
clear 
close all
clc
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
load([pwd,'/Data/UtahArrayData/my147_53'],'LFPs')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess the raw LFP data
flagBandstop = 1 ;
[sigIn,fsTemporal,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

%% bandpass signals
fLow = (15:1:65)' ;
fHigh = (45:1:95)' ;
subBand = [fLow,fHigh] ;
subBand = [30, 80] ;
bandpassSig = find_bandpassSig(sigIn(5,5,:),subBand, fsTemporal, [],0,8,0) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,[],0) ;
% reshapeHilbert = reshape(hilbertSig,100,[]) ;
%absBandpass = abs(reshapeHilbert) ;

%% wavelet transform
freqRange = [30,80] ;
scaleStep = 0.5 ;
wName = 'cmor1.5-1' ;
[cwtCoef,pseudoFreq] = find_cwtCoef2(squeeze(sigIn(5,5,:))', freqRange, scaleStep,...
   fsTemporal,wName) ;
absBandpass = abs(cwtCoef) ;

%% plot the Hilbert spectrum
zscoreCoef = zscore(abs(squeeze( hilbertSig(:,1,1,:) )),1);
figure;
timeStart = 0 ;
timeEnd = 2 ;
plotSection = round(timeStart*fsTemporal+1):round(timeEnd*fsTemporal) ;
timeAxis = timeStart:1/fsTemporal:timeEnd ;
imagesc(timeAxis,30:3:80,(((zscoreCoef(:,plotSection)) )))
% 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(gca,'YDir','normal')

%% plot the wavelet spectrum
timeStart = 29 ;
timeEnd = 31 ;
plotSection = round(timeStart*fsTemporal+1):round(timeEnd*fsTemporal) ;
timeAxis = linspace(timeStart,timeEnd,length(plotSection)) ;

%zscoreCoef = zscore(abs(cwtCoef),1);
zscoreCoef =abs(cwtCoef) ;
figure
% uimagesc(tempTimeAxis,pseudoFreq(2000:-1:20),(abs(wt(2000:-1:20,:)) ))
uimagesc(timeAxis,pseudoFreq(end:-1:1),(((zscoreCoef(end:-1:1,plotSection)) )))
% 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(gca,'YDir','normal')

%% find burst 1D
GammaBurstEvent = find_Burst_1D(absBandpass(end,:),fsTemporal) ;

%% find burst 2D
Y = prctile(zscoreCoef(:),95) ;
binaryImage = zeros(size(zscoreCoef)) ;
binaryImage(zscoreCoef>=Y) = 1 ;
%%
CC = bwconncomp(binaryImage(:,:)) ;
validRegionIdx = 0;
count = 1;
for iRegion = 1:size(CC.PixelIdxList,2)
    if(size(CC.PixelIdxList{iRegion},1)>200)
        validRegionIdx(count) = iRegion ;
        count = count + 1 ;
    end
end
% for iRegion = 1:length(validRegionIdx)
%     figure
%     temp = zeros(size(binaryImage(:,1:2*fsTemporal))) ;
%     temp(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
%     uimagesc(timeAxis(2:end),pseudoFreq(end:-1:1),temp(end:-1:1,:))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
% end
S = regionprops(CC,'Centroid');
centroids = cat(1, S.Centroid);

cTime = centroids(validRegionIdx(:),1)/fsTemporal ;
cFreq = pseudoFreq(round(centroids(validRegionIdx(:),2))) ;

B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

duTime = boundary(validRegionIdx(:),3)/fsTemporal ;
freqLower = min([length(pseudoFreq)*ones(length(validRegionIdx),1),...
    round( centroids(validRegionIdx(:),2)+boundary(validRegionIdx(:),4)/2)]') ;
freqUpper = max([1*ones(length(validRegionIdx),1),...
    round( centroids(validRegionIdx(:),2)-boundary(validRegionIdx(:),4)/2)]') ;
bwFreq = pseudoFreq(freqUpper ) - pseudoFreq(freqLower) ;

%%
figure;
 hist(cFreq,20)
title('histogram of central frequency')
figure;
hist(duTime,40)
title('histogram of duration')
figure;
hist(bwFreq,20)
title('histogram of bandwidth')

figure;
plot(duTime,cFreq,'o')
ylabel('freq(Hz)')
xlabel('duration(s)')

%% find burst 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wavelet transform
timeStart = 1 ;
timeEnd = 11 ;
timeStep = fix(fsTemporal*timeStart+1) : fix(fsTemporal*timeEnd) ;

freqRange = [30,80] ;
scaleStep = 0.5 ;
wName = 'cmor1-5' ;
reshapeSig = reshape(sigIn,100,[]) ;
zscoreCoef = zeros(100,212,length(timeStep)) ;
for ichannel = setdiff(1:size(reshapeSig,1), badChannels)
    [cwtCoef,pseudoFreq] = find_cwtCoef2(reshapeSig(ichannel,timeStep),...
        freqRange, scaleStep,fsTemporal,wName) ;
    zscoreCoef(ichannel,:,:) = zscore(abs(cwtCoef),1);
    disp(['percentage ',int2str(ichannel),'%'])
end
zscoreCoef = reshape(zscoreCoef,10,10,212,length(timeStep)) ;

    %% find burst 2D
Y = prctile(zscoreCoef(:),98) ;
binaryImage = zeros(size(zscoreCoef)) ;
binaryImage(zscoreCoef>Y) = 1 ;
%
CC = bwconncomp(binaryImage) ;
%
validRegionIdx = 0;
count = 1;
for iRegion = 1:size(CC.PixelIdxList,2)
    if(size(CC.PixelIdxList{iRegion},1)>400)
        validRegionIdx(count) = iRegion ;
        count = count + 1 ;
    end
end
%% 
% for iRegion = 1:length(validRegionIdx)
%     figure
%     temp = zeros(size(binaryImage(:,1:2*fsTemporal))) ;
%     temp(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
%     uimagesc(timeAxis(2:end),pseudoFreq(end:-1:1),temp(end:-1:1,:))
%     % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     set(gca,'YDir','normal')
% end
burstSelect = zeros(size(binaryImage)) ;
for iRegion = 1:length(validRegionIdx)
    burstSelect(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
end

[X,Y,Z] = ndgrid(1:size(burstSelect,1),1:size(burstSelect,2),pseudoFreq) ;
figure;pointsize = 30;burstSelect(burstSelect==0) = nan;
temp = burstSelect(:,:,:,2) ;
scatter3(X(:), Y(:), Z(:), pointsize, temp(:));
%%
S = regionprops(CC,'Centroid');
centroids = cat(1, S.Centroid);


% cTime = centroids(validRegionIdx(:),1)/fsTemporal ;
% cFreq = pseudoFreq(round(centroids(validRegionIdx(:),2))) ;
% 
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);
% 
% duTime = boundary(validRegionIdx(:),3)/fsTemporal ;
% freqLower = min([length(pseudoFreq)*ones(length(validRegionIdx),1),...
%     round( centroids(validRegionIdx(:),2)+boundary(validRegionIdx(:),4)/2)]') ;
% freqUpper = max([1*ones(length(validRegionIdx),1),...
%     round( centroids(validRegionIdx(:),2)-boundary(validRegionIdx(:),4)/2)]') ;
% bwFreq = pseudoFreq(freqUpper ) - pseudoFreq(freqLower) ;

%%
cSpace = centroids(validRegionIdx(:),1:2);

plot(cSpace(:,1),cSpace(:,2),'o')
xlabel('ChannelX')
ylabel('ChannelY')
title('Burst Event Distribution over Space')

figure;hist(boundary(validRegionIdx(:),5).*boundary(validRegionIdx(:),6),40)
xlabel('Area(4mmX4mm)')
ylabel('Count')
% figure;
%  hist(cFreq,20)
% title('histogram of central frequency')
% figure;
% hist(duTime,40)
% title('histogram of duration')
% figure;
% hist(bwFreq,20)
% title('histogram of bandwidth')
% 
% figure;
% plot(duTime,cFreq,'o')
% ylabel('freq(Hz)')
% xlabel('duration(s)')

%%
[X,Y,Z] = ndgrid(1:size(burstSelect,1),1:size(burstSelect,2),pseudoFreq(:)) ;
 for iTime = 1:10000
    temp = burstSelect(:,:,:,iTime) ;
scatter3(X(:), Y(:), Z(:), pointsize, temp(:));
title(['time: ',int2str(iTime),'ms'])
pause(0.01)
 end

 %% write 3D video
 randColor = rand(3,length(validRegionIdx) ) ;
 tempSize = 100*length(pseudoFreq) ;
 % Create video file
fig = figure;
vidTitle = 'GammaBurstMoive' ;
vidObj = VideoWriter(vidTitle);
vidObj.FrameRate = 20 ;
open(vidObj);

 for iTime = 1:size(binaryImage,4)
     for region = 1:length(validRegionIdx)       
         tempIdx =CC.PixelIdxList{validRegionIdx(region)}...
             ( (iTime-1)*tempSize<CC.PixelIdxList{validRegionIdx(region)}); 
         tempIdx = tempIdx( tempIdx   < iTime*tempSize) ;
         tempIdx = tempIdx - (iTime-1)*tempSize; 
         if isempty(tempIdx)
             continue
         end
         x = mod(mod(tempIdx,100),10) ;
         y = floor(mod(tempIdx,100)/10) ;
         z = pseudoFreq(floor(tempIdx/100)+1)' ;
%          tri=delaunay(x,y,z);
%          trisurf(tri,x,y,z);
        P = [x,y,z] ;
        k = boundary(P);
        hold on
        trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor',randColor(:,region),'FaceAlpha',0.5,'EdgeAlpha', 0)
         xlim([1 10])
         ylim([1 10])
         zlim([30 80])
         xlabel('ChannelX')
         ylabel('ChannelY')
         zlabel('Freq(Hz)')
         %patch('Faces',tri,'Vertices',[x,y,z])
        view(135,15)
        % axis tight 
     end
     title(['Time: ',num2str(iTime/fsTemporal,'%2.3f'),'s'])
     % pause(0.01)
     % Write to video file
     writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
     cla
     
 end
 % Close video file
close(vidObj);

%% plot in 2D
newImage = zeros(size(binaryImage)) ;
for iRegion = 1:length(validRegionIdx)
    newImage(CC.PixelIdxList{iRegion}) = 1;
end

newReshapeImage = reshape(newImage,100,size(newImage,3),size(newImage,4)) ;
plotImage = zeros(10,10,size(newImage,4)) ;
plotImageTemp = zeros(100,size(newImage,4)) ;

 % Create video file
fig = figure;
vidTitle = 'GammaBurstFullConnect95%' ;
vidObj = VideoWriter(vidTitle);
vidObj.FrameRate = 20 ;
open(vidObj);
for iTime = 1:10:size(binaryImage,4)
    for iFreq = 1% :size(newImage,3)
        ratio = (size(newImage,3)-iFreq)/size(newImage,3)/2+0.5;
        plotImage(:,:,iTime) = ratio*newImage(:,:,iFreq,iTime) ;
        
    end

      for iChannel = setdiff(1:size(reshapeSig,1), badChannels)
          Idx = find(newReshapeImage(iChannel,:,iTime) == 1) ;
          if ~isempty(Idx)
              ratio = (size(newImage,3)-Idx(1))/size(newImage,3)/2+0.5;
            plotImageTemp(iChannel,iTime) = ratio ;
          end
      end
      plotImage = reshape(plotImageTemp,10,10,size(newImage,4)) ;
      
      imagesc(plotImage(:,:,iTime))
      title(['Yellow: 30Hz, Green: 80Hz, Blue: No Burst; Time: ',num2str(iTime/fsTemporal,'%2.3f'), 's'])
      colorbar
      caxis([0,1])
      % pause(0.1)
      writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
      cla
end
% Close video file
close(vidObj);

%% find burst 4D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reshapeZscore = reshape(zscoreCoef,100,size(zscoreCoef,3),size(zscoreCoef,4)) ;
burstImage = zeros(size(reshapeZscore)) ;

for iChannel = setdiff(1:size(reshapeSig,1), badChannels)
    % find burst 2D
    channelData = squeeze(reshapeZscore(iChannel,:,:)) ;
    Y = prctile(channelData,95) ;
    binaryImage = zeros(size(channelData)) ;
    binaryImage(channelData>=Y) = 1 ;
    %
    CC = bwconncomp(binaryImage) ;
    validRegionIdx = 0;
    count = 1;
    newImage = zeros(size(binaryImage)) ;

    for iRegion = 1:size(CC.PixelIdxList,2)
        if(size(CC.PixelIdxList{iRegion},1)>200)
            validRegionIdx(count) = iRegion ;
            count = count + 1 ;
            newImage(CC.PixelIdxList{iRegion}) = 1 ; 
        end
    end
    
    burstImage(iChannel,:,:) = newImage ;
end    
%%
plotImage = zeros(10,10,size(burstImage,3)) ;
count2 = 1;

for iTime = 1:10:size(burstImage,3)
%     for iFreq = 1% :size(newImage,3)
%         ratio = (size(newImage,3)-iFreq)/size(newImage,3)/2+0.5;
%         plotImage(:,:,iTime) = ratio*newImage(:,:,iFreq,iTime) ;
%         
%     end
      plotImageTemp = zeros(100,1) ;

      for iChannel = setdiff(1:size(burstImage,1), badChannels)
          Idx = find(burstImage(iChannel,:,iTime) == 1) ;
          if ~isempty(Idx)
              ratio = (size(burstImage,2)-Idx(1))/size(burstImage,2)/2+0.5;
            plotImageTemp(iChannel) = ratio ;
          end
      end
      plotImage(:,:,count2) = reshape(plotImageTemp,10,10) ;
      imagesc(plotImage(:,:,count2))
      title(['Time: ',num2str(iTime/fsTemporal,'%2.3f'), 's'])
      colorbar
      caxis([0,1])
      pause(0.1)
      cla
      count2 = count2 + 1;
end
plotImage = plotImage(:,:,1:count2) ;
%%
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<5).
opAlpha = 0.2;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields
opBeta = 1;
% Use flag to calculate amplitude velocity fields rather than phase
usePhase = false;

 % Compute optical flow
    disp('Computing optical flow fields...'); tic
    
    % Calculate velocity fields for every trial, and same average number of
    % steps to converge
    % for itrial = 1:size(wvcfs,4)
    [vx, vy, csteps] = opticalFlow2(dist2, badChannels, ...
        opAlpha, opBeta, usePhase);
    vfs = vx + 1i*vy;
    meanCSteps = mean(csteps);
    % fprintf('Processed trial %i\n', itrial)
    % end
    % fprintf('Frequency band %i\n', freqBand)
    toc
    fprintf('Optical flow took %0.1f steps on average to converge.\n', meanCSteps)

    
%% Plot velocity field
for itime = 1:1000
% % phase
% % signalGrid = angle( squeeze(Gamma(:,:,itime)) );
% signalGrid = ( squeeze(Gamma(:,:,itime)) );
% 
% figure
% % Set signal and colors based on whether you are plotting amplitude
% % or phase
% colorMapSpec = pmkmp_new;
% sigLims = [-pi pi];
% vfColor = [0 0 0];
% % Plot signal grid
% imagesc(signalGrid, sigLims)
% colormap(gca, colorMapSpec)

% Amp
% Optionally interpolate spatial grid of data
signalGrid = abs( squeeze(plotImage(:,:,itime)) ) ;

% Set signal and colors based on whether you are plotting amplitude
% or phase
colorMapSpec = parula;
sigLims = [min(abs(signalGrid(:))), max(abs(signalGrid(:)))];
vfColor = [1 1 1];
% Plot signal grid
imagesc(signalGrid, sigLims)
colormap(gca, colorMapSpec)

hold on
%vfScale =  12 ;
%vf = vfs1(:,:,itime,firstBand) * vfScale;
vf = vfs(:,:,itime) * 5;

quiver(linspace(1, size(signalGrid,2), size(vf,2)), ...
    linspace(1, size(signalGrid,1), size(vf,1)), ...
    real(vf), imag(vf), 0, 'Color', vfColor)
set(gca,'YDir','reverse');
hold off

% Update title with current time
title(['Burst Direction at time: ',num2str(itime/fsTemporal,'%2.3f'), 's' ])

pause(0.01)
cla

end