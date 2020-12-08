close all
    plotLength = fix(100*fsTemporal) ;
    numChannel = 40 ;
    sigSmooth = zeros(numChannel,numChannel,plotLength) ;
    % timeStart = 5000 ;
    
    addpath(genpath([pwd,'/ToolNeuroPatt']))
    addpath(genpath([pwd,'/ToolOthers/nanconv']))
    
    resizeScale = 4 ;
    for iTime = 1:plotLength
        timeSlot = iTime+timeStart ;
        %Try smoothing
        filtWidth = 3;
        filtSigma = 0.6;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        smoothTemp = nanconv((squeeze(sigOri(:,:,timeSlot))),imageFilter,'edge', 'nonanout');
        sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);
    end
    clearvars smoothTemp
    
    subBand = [30,80] ;
    bandpassSig = find_bandpassSig(sigSmooth,subBand, fsTemporal,3) ;
    % bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3) ;
    hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;

    sigIn3 = abs(squeeze(hilbertSig(1,:,:,0*fsTemporal+1:85*fsTemporal) ));

%%
timeStart = 60000 ;

close all
figure; imagesc(squeeze(sigIn2(8,:,timeStart:timeStart+500))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_500ms_020.eps'],'epsc')

figure; imagesc(squeeze(sigIn(4,:,timeStart:timeStart+500))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_500ms_010.eps'],'epsc')

figure; imagesc(squeeze(sigIn3(16,:,timeStart:timeStart+500))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_500ms_040.eps'],'epsc')

figure; imagesc(squeeze(sigIn2(8,:,timeStart:timeStart+1000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_1000ms_020.eps'],'epsc')

figure; imagesc(squeeze(sigIn(4,:,timeStart:timeStart+1000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_1000ms_010.eps'],'epsc')

figure; imagesc(squeeze(sigIn3(16,:,timeStart:timeStart+1000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_1000ms_040.eps'],'epsc')

figure; imagesc(squeeze(sigIn2(8,:,timeStart:timeStart+5000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_5000ms_020.eps'],'epsc')

figure; imagesc(squeeze(sigIn(4,:,timeStart:timeStart+5000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_5000ms_010.eps'],'epsc')

figure; imagesc(squeeze(sigIn3(16,:,timeStart:timeStart+5000))')
xlabel('Space (electrode)')
ylabel('Time (ms)')
set(gca,'YDir','normal')
saveas(gcf,[pwd,'/Results/Project1/space_time/TimeStart',int2str(timeStart),'_5000ms_040.eps'],'epsc')
