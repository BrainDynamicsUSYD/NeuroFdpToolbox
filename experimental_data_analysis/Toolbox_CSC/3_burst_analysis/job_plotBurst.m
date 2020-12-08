function job_plotBurst(arrayID)
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
load([pwd,'/Data/UtahArrayData/my147_53'],'LFPs')

% preprocess the raw LFP data
flagBandstop = 1 ;
[sigIn,fsTemporal,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
switch arrayID
case 1
subBand2 = [1, 4];
prcValue = 90 ;
case 2
subBand2 = [0.5, 4];
prcValue = 95 ;
case 3
subBand2 = [0.5, 4];
prcValue = 90 ;
end

subBand = [30 , 100] ;

bandpassSig = find_bandpassSig(sigIn,subBand, fsTemporal) ;
bandpassSig = reshape(bandpassSig,1,10,10,[]) ;
hilbertSig = find_Hilbert(bandpassSig, fsTemporal) ;
amp = abs(squeeze(hilbertSig(1,:,:,11*fsTemporal+1:110*fsTemporal) ));
clearvars bandpassSig hilbertSig
%%
% subBand2 = [1, 4];
bandpassSig2 = find_bandpassSig(sigIn,subBand2, fsTemporal) ;
bandpassSig2 = reshape(bandpassSig2,1,10,10,[]) ;
hilbertSig2 = find_Hilbert(bandpassSig2, fsTemporal) ;

phase = angle(squeeze(hilbertSig2(1,:,:,11*fsTemporal+1:110*fsTemporal) ));

clearvars sigIn bandpassSig2 hilbertSig2
%% Burst movie

close all
plotLength = fix(10*fsTemporal) ;
phaseSmooth = zeros(20,20,plotLength) ;
ampSmooth = zeros(20,20,plotLength) ;
timeStart = 5000 ;

addpath(genpath([pwd,'/ToolOthers/ToolNeuroPatt']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))

resizeScale = 2 ;
for i = 1:plotLength
    timeSlot = i+timeStart ;
     %Try smoothing
    filtWidth = 3;
    filtSigma = 2;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    phaseSmoothTemp = nanconv(abs(squeeze(phase(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    phaseSmooth(:,:,i) = imresize(phaseSmoothTemp, resizeScale);
    
    
             %Try smoothing
    filtWidth = 3;
    filtSigma = 2;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    ampSmoothTemp = nanconv((squeeze(amp(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
    ampSmooth(:,:,i) = imresize(ampSmoothTemp, resizeScale);
   
end
%%
phasePlot = zeros(size(phaseSmooth)) ;
phasePlot(phaseSmooth>pi*2/3) = 1 ;
%phasePlot(phaseSmooth<-pi*11/12) = 1 ;
phasePlot = phasePlot.*phaseSmooth ;

prc95 = prctile(ampSmooth(:),prcValue) ;
ampPlot = zeros(size(ampSmooth)) ;
ampPlot(ampSmooth>prc95) = 1;
ampPlot = ampPlot.*ampSmooth ;
%%
close all
vidTitle = ['/Project1/Results/burstMovies/',num2str(arrayID)] ;
vidObj = VideoWriter(vidTitle);
vidObj.FrameRate = 20 ;
open(vidObj);
fig=figure ;
set(gcf,'Position',[260 23 1159 926])

for i = 1:plotLength
    timeSlot = i+timeStart ;
    subplot(1,2,1)
    colorMapSpec = pmkmp_new;
    sigLims = [-pi pi];
    % Plot signal grid
    imagesc(phasePlot(:,:,i), sigLims)
    set(gca,'YDir','normal')
    title(['1-4Hz Delta phase at ', int2str(timeSlot/1000*1000),'ms'])
    colormap(gca, colorMapSpec)
    colorbar
    
    
    subplot(1,2,2)
    imagesc(ampPlot(:,:,i))
    set(gca,'YDir','normal')
    title(['30-100Hz Gamma amplitude at ', int2str(timeSlot/1000*1000),'ms'])
    caxis([0 max(ampPlot(:))])
    colorbar
    
    writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    %pause
   
end
close(vidObj);
