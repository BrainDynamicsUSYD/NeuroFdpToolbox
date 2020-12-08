timeUnit = 's';
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;
vidTitle = [pwd,'/experiment'] ;
    vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 20 ;
    open(vidObj);
    fig=figure ;
    % set(gcf,'Visible','On','Position',[104 443 1880 455])
    set(gcf,'Visible','On','Position',[31 299 1525 587])
    set(gcf,'Visible','On','Position',[360 30 790 665])
    
%for iTime = 1:10:20000;
DeltaMin = min(min(min(abs(hilbertSig(1,:,:,:))))) ;
DeltaMax = 0.4*max(max(max(abs(hilbertSig(1,:,:,:)))))  ;

ThetaMin = min(min(min(abs(hilbertSig(2,:,:,:))))) ;
ThetaMax = 0.4*max(max(max(abs(hilbertSig(2,:,:,:)))))  ;

GammaMin = min(min(min(abs(hilbertSig(5,:,:,:))))) ;
GammaMax = 0.4*max(max(max(abs(hilbertSig(5,:,:,:)))))  ;

oriMin =  min(abs(hilbertOri(:))) ;
oriMax =  0.4*max(abs(hilbertOri(:)))  ;

timeDiff = fix(1*fsTemporal) ;

for timeHilbert = 1:1:2000;
    % timeHilbert = floor(iTime/10) +1 ;
    subplot(2,2,1);
    imagesc(abs(squeeze(hilbertSig(1,:,:,timeHilbert)))); 
    caxis([DeltaMin,DeltaMax]); 
        title(['Delta signals (1-4Hz) at',num2str(timeHilbert),'ms']);

    subplot(2,2,2);imagesc(abs(squeeze(hilbertSig(2,:,:,timeHilbert))));
    caxis([ThetaMin,ThetaMax ]);
        title(['Theta signals (5-8Hz) at',num2str(timeHilbert),'ms']);

    subplot(2,2,3);imagesc(abs(squeeze(hilbertSig(5,:,:,timeHilbert))));
    caxis([GammaMin,GammaMax ]);
        title(['Gamma signals (30-80Hz) at',num2str(timeHilbert),'ms']);

%     subplot(2,2,4);imagesc(tempI2(:,:,11600+iTime));
%     caxis([min(min(min(abs(tempI2(:,:,:))))) ...
%         0.8*max(max(max(abs(tempI2(:,:,:))))) ]);
%     title(['original signals (sum of abs current) at',num2str(timeHilbert),'ms']);
    subplot(2,2,4);imagesc(squeeze(abs(hilbertOri(1,:,:,timeDiff+timeHilbert))));
    caxis([oriMin, oriMax ]);
    title(['original signals at',num2str(timeHilbert),'ms']);
    
    writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    % pause
    cla;
end

close(vidObj);
