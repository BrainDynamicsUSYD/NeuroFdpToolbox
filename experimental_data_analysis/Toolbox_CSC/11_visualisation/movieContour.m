sigIn = abs(sigSmooth) ;
sigMin = min(sigIn(:)) ;
sigMax = max(sigIn(:)) ;
prc95 =    prctile(sigIn(:), 95) ;
prc99 =    prctile(sigIn(:), 99)  ;




for iTime= 200: 1000;
    c{iTime} = contourc(sigIn(:,:,iTime),[prc95,0]) ;
    
end


% 
% for iTime= 200: 1000;
%     c1{iTime} = contourc(sigIn(:,:,iTime),[prc95]) ;
%     
% end
% for iTime= 200: 1000;
%     c2{iTime} = contourc(sigIn(:,:,iTime),[prc99]) ;
%     
% end
% for iTime= 200: 1000;
%     plot(c1{iTime})
%     pause
%     cla
% end
%%
close all
% timeUnit = 's';
% t = datetime('now') ;
% dateStr = datestr(t,'mmmmdd_HH:MM') ;
% % vidTitle = [pwd,'/experiment'] ;
% vidTitle = ['BetaMovie_',dateStr ];
%     vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
%     v.Quality = 50 ;
%     vidObj.FrameRate = 20 ;
%     open(vidObj);

timeStart = 50000 ;
timeEnd = 60000 ;
for iTime= timeStart: timeEnd;
    c{iTime} = contourc(sigIn(:,:,iTime),[prc95,0]) ;
    
end

    fig=figure ;
    % set(gcf,'Visible','On','Position',[104 443 1880 455])
    set(gcf,'Visible','On','Position',[31 299 1525 587])
    set(gcf,'Visible','On','Position',[676 357 834 746])
    
for iTime= timeStart: timeEnd;
    imagesc(1:0.5:10,1:0.5:10,sigIn(:,:,iTime))
    colorbar
    caxis([sigMin 0.3*sigMax])
    title(['Beta signals (13-30Hz) at ',num2str(iTime),'ms (black: 95 prctile)']);
    hold on;
    plot((c{iTime}(1,2:end)+0.5)/2,(c{iTime}(2,2:end)+0.5)/2,'k.','markersize',12)
    hold on;
    % plot(c2{iTime}(1,2:end)/2-0.5,c2{iTime}(2,2:end)/2-0.5,'r.','markersize',12)
    xlim([1,10])
    ylim([1,10])
pause
    % writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    % pause
    cla;
end

% close(vidObj);
