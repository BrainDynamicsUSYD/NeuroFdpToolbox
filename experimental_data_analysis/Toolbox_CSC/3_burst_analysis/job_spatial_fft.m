function job_spatial_fft(arrayID)
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
cd ..
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
%%
arrayID = 8 ;
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
        subBand = [30 , 100] ;0.5
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
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
clearvars LFPs
sigOriTemp = reshape(sigOri,10*10,[]) ;

%%
sigOriComplex = nan(100,size(sigOriTemp,2)) ;
for iChannel = setdiff(1:100,badChannels)
    sigOriComplex(iChannel,:) = hilbert(sigOriTemp(iChannel,:)) ;
end
sigHilbert = reshape(sigOriComplex,10,10,[]) ;

%%
% vidTitle = [pwd,'/Results/RawSignalMovies/',dataFileName,'_',...
%         '.avi'] ;
%     vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);
     fig=figure ;
     set(gcf,'Position',[260 23 1159 926])     
    plotAmp = 0 ;
    for iTime = 1:1000
        
        imagesc(angle(sigHilbert(:,:,iTime)))
        set(gca,'YDir','normal')
        title(['Raw signal at ', ...
            int2str(iTime/fsTemporal*1000),'ms'])
        if plotAmp
            caxis([0 0.6*max(abs(sigHilbert(:)))])
        else
            colorMapSpec = pmkmp_new;
            sigLims = [-pi pi];
            colormap(gca, colorMapSpec)
            caxis(sigLims)
        end
        colorbar
        
   
        % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
        pause(0.01)
        cla
        
    end
    close(vidObj);
    
%% Interpolation and smoothing

close all
% plotLength = fix(100*fsTemporal) ;
% sigSmooth = zeros(20,20,size(sigHilbert,3)) ;
% sigSmooth2 = zeros(20,20,size(sigHilbert,3)) ;
% % timeStart = 5000 ;
% 
addpath(genpath([pwd,'/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))
% 
% resizeScale = 2 ;
for iTime = 1:10000 %size(sigIn,3)% plotLength
    timeSlot = iTime ;
     %Try smoothing
    filtWidth = 5;
    filtSigma = 0.2;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp(:,:,timeSlot) = nanconv(abs(squeeze(sigHilbert(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    % smoothTemp = imgaussfilt(abs(squeeze(sigHilbert(:,:,timeSlot))), filtSigma ) ;
    
    filtWidth = 5;
    filtSigma = 0.4;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp2(:,:,timeSlot) = nanconv(abs(squeeze(sigHilbert(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    % smoothTemp2 = imgaussfilt(abs(squeeze(sigHilbert(:,:,timeSlot))), filtSigma ) ;
    
    filtWidth = 5;
    filtSigma = 0.8;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp3(:,:,timeSlot) = nanconv(abs(squeeze(sigHilbert(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    % smoothTemp3 = imgaussfilt(abs(squeeze(sigHilbert(:,:,timeSlot))), filtSigma ) ;
    
    filtWidth = 5;
    filtSigma = 1.6;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp4(:,:,timeSlot) = nanconv(abs(squeeze(sigHilbert(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    % smoothTemp4 = imgaussfilt(abs(squeeze(sigHilbert(:,:,timeSlot))), filtSigma ) ;
    
    filtWidth = 5;
    filtSigma = 3.2;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp5(:,:,timeSlot) = nanconv(abs(squeeze(sigHilbert(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    % smoothTemp5 = imgaussfilt(abs(squeeze(sigHilbert(:,:,timeSlot))), filtSigma ) ;
end

%%
timeSelect = 1 ;
figure;
imagesc(squeeze(abs(sigHilbert(:,:,timeSelect))))
figure;
imagesc(squeeze(smoothTemp(:,:,timeSelect)))
figure;
imagesc(squeeze(smoothTemp2(:,:,timeSelect)))
figure;
imagesc(squeeze(smoothTemp3(:,:,timeSelect)))
figure;
imagesc(squeeze(smoothTemp4(:,:,timeSelect)))
figure;
imagesc(squeeze(smoothTemp5(:,:,timeSelect)))

%%
Bandpass1 = smoothTemp - smoothTemp2 ;
Bandpass2 = smoothTemp2 - smoothTemp3 ;
Bandpass3 = smoothTemp3 - smoothTemp4 ;
Bandpass4 = smoothTemp4 - smoothTemp5 ;

%%
figure;
imagesc(squeeze(abs(sigHilbert(:,:,timeSelect))))
figure;
imagesc(squeeze(Bandpass1(:,:,timeSelect)))
figure;
imagesc(squeeze(Bandpass1(:,:,timeSelect)))
figure;
imagesc(squeeze(Bandpass3(:,:,timeSelect)))
figure;
imagesc(squeeze(Bandpass4(:,:,timeSelect)))


%%
P = 95 ;
S1 = prctile(Bandpass1(:),P) ;
S2 = prctile(Bandpass2(:),P) ;
S3 = prctile(Bandpass3(:),P) ;
S4 = prctile(Bandpass4(:),P) ;

%% 
eddie1 = zeros(size(Bandpass1)) ;
eddie1(Bandpass1>=S1) = 1 ;

eddie2 = zeros(size(Bandpass2)) ;
eddie2(Bandpass2>=S2) = 1 ;

eddie3 = zeros(size(Bandpass3)) ;
eddie3(Bandpass3>=S3) = 1 ;

eddie4 = zeros(size(Bandpass4)) ;
eddie4(Bandpass4>=S4) = 1 ;

%%

fig = figure;
set(gcf,'Position',[260 23 1159 926])     

for iTime = 1:size(eddie1,3)
    subplot(2,2,1)
    imagesc(eddie1(:,:,iTime))
    
    subplot(2,2,2)
    imagesc(eddie2(:,:,iTime))
    
    subplot(2,2,3)
    imagesc(eddie3(:,:,iTime))
    
    subplot(2,2,4)
    imagesc(eddie4(:,:,iTime))
    
    pause
    cla
end


