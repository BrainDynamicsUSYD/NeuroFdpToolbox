function job_burstTrajectory(arrayID)
%% Main Function for spatial temporal patterns
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
%
% Initialisation
% clear 
% close all
% clc
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))

%%

switch arrayID
    case 1
        dataFileName = 'my144_101' ;
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 2
        dataFileName = 'my144_101' ;
        timeStart2 = 5000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;
        
        makeMovie = 1 ;
        
    case 3
        dataFileName = 'my144_101' ;
        timeStart2 = 9000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 4
        dataFileName = 'my147_02' ;
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 5
        dataFileName = 'my147_02' ;
        timeStart2 = 5000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;
        
        makeMovie = 1 ;
        
    case 6
        dataFileName = 'my147_02' ;
        timeStart2 = 9000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
    case 7
        dataFileName = 'ma027_032' ;
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 8
        dataFileName = 'ma027_032' ;
        timeStart2 = 5000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;
        
        makeMovie = 1 ;
        
    case 9
        dataFileName = 'ma027_032' ;
        timeStart2 = 9000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 10
        dataFileName = 'ma025_03' ;
        timeStart2 = 0;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
        
    case 11
        dataFileName = 'ma025_03' ;
        timeStart2 = 5000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;
        
        makeMovie = 1 ;
        
    case 12
        dataFileName = 'ma025_03' ;
        timeStart2 = 9000;
        subBand = [30 , 100] ;
        plotAmp = 1 ;
        plotTrajectory = 0 ;
        % timeStep = 5 ;

        makeMovie = 1 ;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')
%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigIn,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
clearvars LFPs

% find subband signals
bandpassSig = find_bandpassSig(sigIn,subBand, fsTemporal) ;
bandpassSig = reshape(bandpassSig,1,10,10,[]) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal) ;

if plotAmp
    sigIn = abs(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:110*fsTemporal) ));
else
    sigIn = angle(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:110*fsTemporal) ));
end

clearvars bandpassSig hilbertSig


%% Interpolation and smoothing

close all
plotLength = fix(10*fsTemporal) ;
sigSmooth = zeros(20,20,plotLength) ;
timeStart = 5000 ;

addpath(genpath([pwd,'/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))

resizeScale = 2 ;
for iTime = 1:plotLength
    timeSlot = iTime+timeStart ;
     %Try smoothing
    filtWidth = 3;
    filtSigma = 1;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);   
end
clearvars smoothTemp
%% find 95 percentile
if ~plotAmp
    sigBinary = zeros(size(sigSmooth)) ;
    sigBinary(sigSmooth>pi*5/6) = 1 ;
    sigPlot = sigBinary.*sigSmooth ;
else
    prc95 = prctile(sigSmooth(:),95) ;
    sigBinary = zeros(size(sigSmooth)) ;
    sigBinary(sigSmooth>prc95) = 1;
    sigPlot = sigBinary.*sigSmooth ;
end

%% find central
largestCentroids = nan(size(sigBinary,3),2) ;
% count = 1;
for iTime = 1:size(sigBinary,3)
    CC = bwconncomp(sigBinary(:,:,iTime)) ;

    if (size(CC.PixelIdxList,2) == 0)
        continue
    else
        sizeRegion = zeros(size(CC.PixelIdxList,2),1) ;
        for iRegion = 1:size(CC.PixelIdxList,2)
            sizeRegion(iRegion) = size(CC.PixelIdxList{iRegion},1) ;
        end
        [regionSize,largestRegionIdx] = max(sizeRegion) ;
        if regionSize<4
            continue
        else
            S = regionprops(CC,'Centroid');
            Centroids = cat(1, S.Centroid);
            
            largestCentroids(iTime,:) = Centroids(largestRegionIdx,:) ;
            % count = count+1 ;
        end
    end
end

%% calculate the duration of each cluster
plot(largestCentroids(:,1))
hold on
% plot(largestCentroids(:,2))-20
% first find drift location
diffCent = diff(largestCentroids,1) ;
driftLoc = find((abs(diffCent(:,1))>10) | (abs(diffCent(:,2))>10) )  ;

% find continuous bursts
binaryCentTemp = zeros(size(largestCentroids,1),1) ;
binaryCentTemp(~isnan(largestCentroids(:,1))) = 1 ;
binaryCent = binaryCentTemp ;

for iDrift = 1:length(driftLoc)
    binaryCent = [binaryCent(1:(driftLoc(iDrift)+iDrift-1));0;...
        binaryCent((driftLoc(iDrift)+iDrift):end)] ;
end
plot(binaryCent)

%%
CC = bwconncomp(binaryCent) ;
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);
hist(boundary(:,4),40)

%% Centroids with long duration
% first find drift location
diffCent = diff(largestCentroids,1) ;
driftLoc = find((abs(diffCent(:,1))>4) | (abs(diffCent(:,2))>4) )  ;
% find continuous bursts
binaryCentTemp2 = zeros(size(largestCentroids,1),1) ;
binaryCentTemp2(~isnan(largestCentroids(:,1))) = 1 ;
binaryCentTemp2(driftLoc) = 0 ;

CC = bwconncomp(binaryCentTemp2) ;
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

minDuration = 30/1000*fsTemporal ;
selectCentroids = nan(size(largestCentroids)) ;
validIdx = find(boundary(:,4)>minDuration) ;
for idx = 1:length(validIdx)
    selectCentroids(CC.PixelIdxList{validIdx(idx)},:) = 1 ;
end
largestCentroids = largestCentroids.*selectCentroids ;


%% Calculate the distance of Gamma burst
% electrode distance 0.4mm, 0.2mm after interpolation 
distCent = 0.5*sqrt(sum((largestCentroids(2:end,:) - largestCentroids(1:end-1,:)).^2,2)) ;
% distCent = (largestCentroids(2:end,2) - largestCentroids(1:end-1,2)) ;
distCentMM = 0.4*distCent ;
% hist(distanceCent,40)
% xlabel('Distance (mm)')
%% Histogram of Gamma drift distsance
close all
[n, xout] = hist(distCentMM,200);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% xlabel('Distance (unit electrode)')
xlabel('Distance (mm)')
%xlim([0.1,15])
xlim([0,14*0.4])
hold on

% % fit a levy distribution
% pd = fitdist(distCent,'stable') ;
% x_values = -10:1:10;
% y = pdf(pd,x_values);
% plot(x_values,y,'r','LineWidth',2)

%% plot the drift trajectory
if plotAmp
    fileName = [dataFileName,'amp',num2str(subBand(1)),'_', num2str(subBand(2)),'Hz' ] ;
else
    fileName = [dataFileName,'phase',num2str(subBand(1)),'_',num2str(subBand(2)), 'Hz' ] ;
end

if plotTrajectory
    close all
    fig = figure;
    for timeSlot = 0:timeStep:7*timeStep ;
        subplot(4,2,timeSlot/timeStep+1)
        timeScale = fix(timeSlot*fsTemporal)+1 : fix((timeSlot+timeStep)*fsTemporal) ;
        for iTime = timeScale
            currentPoint = largestCentroids(iTime,:) ;
            nextPoint = largestCentroids(iTime+1,:) ;
            plot(currentPoint(:,1),currentPoint(:,2),'.','MarkerSize',20)
            xlim([0,20])
            ylim([0,20])
            hold on
            trajectoryX = linspace(currentPoint(1),nextPoint(1),100) ;
            trajectoryY = linspace(currentPoint(2),nextPoint(2),100) ;
            plot(trajectoryX,trajectoryY)
            xlim([0,20])
            ylim([0,20])
            %pause(0.1)
        end
        title(['time start at:',num2str(timeSlot),' s'])
        hold off
    end
    set(gcf,'Position',[625 35 831 926]) ;
    %print(fig,[pwd,'/Results/GammaBurst/Trajectory/',fileName,...
    %    'Trajectory_',num2str(timeStep),'sStep'],'-dpng')
end
%% movie

%% make a movie
if makeMovie
    vidTitle = [pwd,'/Results/GammaBurst/','Full',fileName,num2str(timeStart2)] ;
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = 20 ;
    open(vidObj);
    fig=figure ;
    set(gcf,'Position',[260 23 1159 926])
    
    
    for iTime = 1:1000
        timeSlot = iTime+timeStart2 ;
        subplot(2,2,1)
        imagesc(sigIn(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        title(['Raw signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,2)
        imagesc(sigSmooth(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Smoothed signal at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,3)
        imagesc(sigPlot(:,:,iTime))
        set(gca,'YDir','normal')
        title(['Detected burst at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/1000*1000),'ms'])
        if plotAmp
            caxis([0 0.5*max(sigIn(:))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
        subplot(2,2,4)
        imagesc(sigSmooth(:,:,iTime))
        currentPoint = largestCentroids(iTime,:) ;
        % nextPoint = largestCentroids(iTime+1,:) ;
        plot(currentPoint(:,1),currentPoint(:,2),'.','MarkerSize',20)
        xlim([0,20])
        ylim([0,20])
        title(['Detected burst Centre at ',num2str(subBand),'Hz at ', ...
            int2str(timeSlot/fsTemporal*1000),'ms'])
        colorbar
        
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        pause(0.01)
        
    end
    close(vidObj);
end

% %% 5 by 2 smooth images
% close all
% figure;
% stepSize = 5 ;
% count = 1 ; 
% timeStart = 1020 ;
% for iTime = timeStart:stepSize:(stepSize*9+timeStart)
%     subplot(2,5,count)
%     count = count+1 ;
%     
%     imagesc(sigSmooth(:,:,iTime))
%     set(gca,'YDir','normal')
%     caxis([0,0.3*max(sigSmooth(:))])
% end