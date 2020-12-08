function moviesP2

%% movies
flagMovie1 = 0 ;
flagMovie2 = 0 ;
flagMovie3 = 0 ;


%% movie of different bands
if flagMovie1
    for iTime = 1000:2000
        iImage = 1 ;
        for iBand = 11:2:19
        sigPlot = sigIn{iBand}(:,:,iTime) ;
        subplot(2,3,iImage)
        imagesc(sigPlot)
        iImage = iImage + 1 ;
        end
        pause(0.1)
        cla
    end

end


%% movie of spatial temporal shift patterns
if flagMovie2
    
    for iNum =1:numBand
        maxSize = [] ;
        
        fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
        if surFlag == 0
        load([saveFolderName,fileName(1).name],'burstDetect')
        else
            load([saveFolderName,fileName(2).name],'burstDetect')
        end
        binaryBurst = zeros(size(sigIn{1})) ;   % or wvcfs
        for iBurst = 1:length(burstDetect.PixelIdxList)
            binaryBurst(burstDetect.PixelIdxList{iBurst}) = 1 ;
        end    
        sigPlot{iNum} = sigIn{iNum}.*binaryBurst ;
    end

deltaT = [0,64,4096] ;
deltaX = [0,4,16] ;
deltaY = [0,4,16] ;
sizeDelS = length(deltaX)*length(deltaY) ;
sizeDelT = length(deltaT) ;
deltaX = meshgrid(deltaX) ;
deltaY = meshgrid(deltaY)' ;

coef3 = zeros(numBand,numBand,sizeDelS*sizeDelT) ;

tic
timeLength = size(sigIn{1},3) ;
for iDelT = 1:sizeDelT
    for iDelS = 1:sizeDelS
        for iBand1 = [4,13] %1:numBand
            for iBand2 = [4,13] %1:numBand
                %[227 297 1437 558]
                
                b1_x = 1:spaceLength-deltaX(iDelS) ;
                b1_y = 1:spaceLength-deltaY(iDelS) ;
                %b1 = bsxfun(@minus,sigIn{iBand1}(b1_x,b1_y,1:timeLength-deltaT(iDelT)) , mean(mean(...
                %    sigIn{iBand1}(:,:,1:timeLength-deltaT(iDelT))))) ;
                %             b1 = (wvcfs{iBand1}(:,:,1:timeLength-deltaT(iDelT))) - mean(mean((...
                %                 wvcfs{iBand1}(:,:,1:timeLength-deltaT(iDelT))))) ;
                image1 = sigIn{iBand1}(b1_x,b1_y,1:timeLength-deltaT(iDelT)) ;
                
                b2_x = deltaX(iDelS)+1:spaceLength ;
                b2_y = deltaY(iDelS)+1:spaceLength ;
                %b2 = bsxfun(@minus,sigIn{iBand2}(b2_x,b2_y,deltaT(iDelT)+1:timeLength) , mean(mean(...
                %    sigIn{iBand2}(:,:,deltaT(iDelT)+1:timeLength)))) ;
                %             b2 = (wvcfs{iBand2}(:,:,deltaT(iDelT)+1:timeLength)) - mean(mean((...
                %                 wvcfs{iBand2}(:,:,deltaT(iDelT)+1:timeLength)))) ;
                image2 = sigIn{iBand2}(b2_x,b2_y,deltaT(iDelT)+1:timeLength) ;
                
                subplot(1,2,1)
                imagesc(image1)
                subplot(1,2,2)
                imagesc(image2)
                pause
                cla
                %b1b1 = bsxfun(@times,b1,b1) ;
                %b2b2 = bsxfun(@times,b2,b2) ;
                %b1b2 = bsxfun(@times,b1,b2) ;
                %coef3(iBand1,iBand2,iDelS+(iDelT-1)*sizeDelS) = sum((sum(sum(b1b2))))/ sqrt...
                %    (sum(((sum(sum(b1b1))))*sum((sum(sum(b2b2)))))) ;
            end
        end
        disp(['Space shift finish ',num2str(iDelS/sizeDelS*100),' %' ]);
    end
    disp(['Time delay finish ',num2str(iDelT/sizeDelT*100),' %' ]);
    toc
end

end


%% 3D representation of patterns
if flagMovie3
close all
    for iNum =5:4:17
        maxSize = [] ;
        
        fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
        if surFlag == 0
        load([saveFolderName,fileName(1).name],'burstDetect')
        else
            load([saveFolderName,fileName(2).name],'burstDetect')
        end
        binaryBurst = zeros(size(sigIn{1})) ;   % or wvcfs
        for iBurst = 1:length(burstDetect.PixelIdxList)
            binaryBurst(burstDetect.PixelIdxList{iBurst}) = 1 ;
        end    
        sigPlot{iNum} = sigIn{iNum}.*binaryBurst ;
    end
    countBand = 1 ;
    for iNum =5:4:17
        sigTemp{countBand} = sigPlot{iNum}(:,:,8800:9200) ;  % 8000:12000
        countBand = countBand+1 ;
    end
    for iNum2 = 1:4 
        figure;
        set(gcf,'Position',[435*(iNum2-1) 102 435 640])
        [x,y,z] = meshgrid(1:20,1:20,1:size(sigTemp{iNum2},3)) ;
        data = zscore(sigTemp{iNum2}) ;
        p = patch(isosurface(data,mean(data(:))));
        isonormals(data,p)
        [r,g,b] = meshgrid(20:-1:1,1:20,1:size(sigTemp{iNum2},3));
        c = isocolors(r/20,g/20,b/size(sigTemp{iNum2},3),p);
        p.FaceVertexCData = 1-c;
        p.FaceColor = 'interp';
        p.EdgeColor = 'none';
        xlabel('xAxis (electrodes)')
        ylabel('yAxis (electrodes)')
        zlabel('time (ms)')
        % view(150,30)
        view(-14,3)
        % daspect([1 1 1])
        camlight
        lighting gouraud
    end


end

%% 2D representation
if flagMovie4
close all
    %for iNum =5:4:17
    iNum = 17 ;
        maxSize = [] ;
        
        fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
        if surFlag == 0
        load([saveFolderName,fileName(1).name],'burstDetect')
        else
            load([saveFolderName,fileName(2).name],'burstDetect')
        end
        binaryBurst = zeros(size(sigIn{1})) ;   % or wvcfs
        for iBurst = 1:length(burstDetect.PixelIdxList)
            binaryBurst(burstDetect.PixelIdxList{iBurst}) = 1 ;
        end    
        sigPlot{iNum} = sigIn{iNum}.*binaryBurst ;
%         WCentroids = zeros(size(binaryBurst,3),2) ;
%         for iTime = 1:size(binaryBurst,3)
%             S = regionprops(binaryBurst(:,:,iTime),sigPlot{iNum}(:,:,iTime),{'Centroid','WeightedCentroid'} );
%             if isempty(cat(1, S.WeightedCentroid))
%                 continue
%             end
%             WCentroids(iTime,:) = cat(1, S.WeightedCentroid) ; 
%         end
 %%          
        for iBurst = 1:length(burstDetect.PixelIdxList)  %2Hz: 41, 57
           timeStart = fix(burstDetect.PixelIdxList{iBurst}(1)/20/20)+fix(1/cfreqso(iNum)*1000/2) ;
           subplot(1,2,1)
           imagesc(sigPlot{iNum}(:,:,timeStart) )
%           hold on
%           plot(WCentroids(timeStart,1),WCentroids(timeStart,2),'r.','markersize',8)
           
           subplot(1,2,2)
           imagesc(sigPlot{iNum}(:,:,timeStart+fix(1/cfreqso(iNum)*1000)) ) %1cycle
%           hold on
%            for iTime = timeStart:timeStart+fix(1/cfreqso(iNum)*1000)
%             plot([WCentroids(iTime,1),WCentroids(iTime+1,1)],[WCentroids(iTime,2),...
%                 WCentroids(iTime+1,2)],'r.-','markersize',8)
%            end
           set(gcf,'Position',[675 625 877 336])
           pause
           cla
        end
    %end
    %%
    
    for iTime = 1:10000
            countBand = 1;
        sigTemp{countBand} = sigPlot{iNum}(:,:,iTime) ;  % 8000:12000
        subplot(4,1,countBand)
        imagesc(sigTemp{countBand})
        countBand = countBand+1 ;
        
    pause
    cla
    end
      
end
%% spatial shifts and overlap
saveFolderName = [pwd,'/Results_data/Project2/BurstDetection6/'] ;

iNum = 5 ; 
fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
if surFlag == 0
    load([saveFolderName,fileName(1).name],'burstDetect')
else
    load([saveFolderName,fileName(2).name],'burstDetect')
end
% sizeData = size(sigIn{1}) ;
sizeData = [20,20,4e5] ;
binaryBurst = zeros(sizeData) ;   % or wvcfs
for iBurst = 1:length(burstDetect.PixelIdxList)
    binaryBurst(burstDetect.PixelIdxList{iBurst}) = 1 ;
end

iNum = 13 ; 
fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
if surFlag == 0
    load([saveFolderName,fileName(1).name],'burstDetect')
else
    load([saveFolderName,fileName(2).name],'burstDetect')
end
binaryBurst2 = zeros(sizeData) ;   % or wvcfs
for iBurst = 1:length(burstDetect.PixelIdxList)
    binaryBurst2(burstDetect.PixelIdxList{iBurst}) = 1 ;
end

%%
for deltaX = 1:20;
% deltaX = 16;
for deltaY = 1:20 ;
% deltaY = 8 ;

spaceLength = 20 ;
b1_x = 1:spaceLength-deltaX ;
b1_y = 1:spaceLength-deltaY ;
binaryBurstShift = binaryBurst(b1_x,b1_y,:) ;

b2_x = deltaX+1:spaceLength ;
b2_y = deltaY+1:spaceLength ;
binaryBurstShift2 = binaryBurst2(b2_x,b2_y,:) ;

CC_burstShift = bwconncomp(binaryBurstShift) ;    
CC_burstShift2 = bwconncomp(binaryBurstShift2) ;    

%% overlap
% 1 and its overlap
R = [] ;
for iBurst = 1:CC_burstShift.NumObjects
    R(iBurst) = length(find(binaryBurstShift2(CC_burstShift.PixelIdxList{iBurst})==1))...
        /length(CC_burstShift.PixelIdxList{iBurst}) ;
end
R_mean(deltaX,deltaY) = mean(R) ;
% 2 and its overlap
R2 = [] ;
for iBurst = 1:CC_burstShift2.NumObjects
    R2(iBurst) = length(find(binaryBurstShift(CC_burstShift2.PixelIdxList{iBurst})==1))...
        /length(CC_burstShift2.PixelIdxList{iBurst}) ;
end
R_mean2(deltaX,deltaY) = mean(R2) ;
end
end