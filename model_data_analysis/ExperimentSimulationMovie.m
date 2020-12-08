% make a movie
close all
% sigBurst = sigIn.*burstIdx ;
makeMovie = 1;
if makeMovie
%     vidTitle = [pwd,'/Results/movies/New',dataFileName,'_', num2str(boundPrctile) ,'%_',...
%         num2str(minBurstTime),'ms',num2str(timeStart2),'.avi'] ;
%     vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);
vidObj = VideoWriter('ExperimentSimulation.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
    fig = figure;
     set(gcf,'Position',[294 459 1169 453])   
     period = 7230:7280;
    temp = LFPr(:,:,period+1.2*fs);
    mm = minmax(temp(:)');

    for iTime = 1700:1750% 1600:2000
        subplot(1,2,1)

        imagesc(sigSmooth(:,:,iTime))       
        % imagesc(sigBurst(:,:,iTime+timeStart))
        set(gca,'YDir','normal')
        %title(['Detected burst at ',num2str(subBand),'Hz at ', ...
        %    int2str(timeSlot/fsTemporal*1000),'ms'])
            caxis([0 0.35*max(sigSmooth(:))])

        colorbar
        ts = sprintf('t = %8.1f ms',iTime);
        title(ts);
        

   subplot(1,2,2)
        LFP = LFPr(:,:,iTime+5530+1.2*fs);
        imagesc(flipud(LFP),mm)
         colorbar
        ts = sprintf('t = %8.1f ms',iTime);
        title(ts);
        
        
        F = getframe(fig);
            writeVideo(vidObj,F.cdata);
        pause(0.01)
%         cla
        
    end
    close(gcf);
close(vidObj);
    % close(vidObj);
end