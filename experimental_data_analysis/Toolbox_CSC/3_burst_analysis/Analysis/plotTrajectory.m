for j = 38: size(WCentroids,2)
center = WCentroids{j} ;
%for i = 1:size(center,1)
plot(center(:,1),center(:,2),'.-','markersize',20); hold on
%end
xlim([0,10])
ylim([0,10])
title(num2str(j))
pause
close all
end

%%
for i = 1:8 ;  
    subplot(2,4,i);  
    imagesc(1:10,1:10,sigSmoothBackup(:,:,7065+i*20)); 
    set(gca,'YDir','normal'); 
    caxis([0 0.3*max(sigSmoothBackup(:))]) 
end
