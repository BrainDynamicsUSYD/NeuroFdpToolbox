% This is the main function of the project Gamma oscillations organize as 
% transient super-diffusive patterns (experimental part)
%
% Author: Xian Long,  Supervisor: Pulin Gong
%
% Date: 30/04/2019
%
function Project1(arrayID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation
% clear ; 
close ;
clc ;
cd ..
%% add path
addpath(genpath([pwd,'/Toolbox_CSC']))         % Tool
addpath(genpath([pwd,'/ToolOthers/uimage']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))
addpath(genpath([pwd,'/ToolOthers/UTILITIES_Paul']))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch arrayID
    case 1
        surMethodNum = 1 ;
    case 2
        surMethodNum = 2 ;
    case 3
        surMethodNum = 3 ;
    case 4
        surMethodNum = 4 ;
    case 5
        surMethodNum = 5 ;
    case 6
        surMethodNum = 6 ;
    case 7
        surMethodNum = 7 ;
    case 8
        surMethodNum = 8 ;
    case 9
        surMethodNum = 9 ;
    case 10
        surMethodNum = 10 ;
    case 11
        surMethodNum = 0 ;
end
flagPreAnalysis = 0 ;
flagPostAnalysis = 0 ;
%% Load in data
dataFileName = 'ma027_032' ;
load([pwd,'/Data/UtahArrayData/',dataFileName],'LFPs','Fs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing (optional, includes re-organizing and bandstop electric lines)
fsTemporal = Fs ;
flagBandstop = 1 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate surrogate data (optional)
surSig = generateSur(sigOri,surMethodNum,badChannels) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find subband signals
% apply bandpass filtering for Gamma signals
subBand = [30,80] ;
bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels) ;
% Hilbert transform for analytic signals
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
% find amplitdue of the analytic Gamma as the input
sigIn = abs(squeeze((hilbertSig))) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial rescaling and dead electrode interpolation (optional)
% for fig.1E and fig.2
resizeScale = 8 ;
% sigSmooth = sigIn ;

% resizeScale = 2 ;
numChannelS = resizeScale * size(sigIn,1);
sigSmooth = zeros(numChannelS,numChannelS,size(sigIn,3)) ;

for iTime = 1:size(sigIn,3)
    %Try smoothing
    filtWidth = 3;
    filtSigma = 0.6;   % 0.6
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    smoothTemp = nanconv(abs(squeeze(sigIn(:,:,iTime))),imageFilter,'edge', 'nonanout');
    sigSmooth(:,:,iTime) = imresize(smoothTemp, resizeScale);
end
clearvars smoothTemp
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdVal = 2.5;

if flagPreAnalysis
    preAnalysis(sigOri,bandpassSig,sigIn,sigSmooth,badChannels,fsTemporal,stdVal) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst patterns detection
numChannelS = size(sigSmooth,1)*size(sigSmooth,2) ;
sigReshape = reshape(sigSmooth,numChannelS,[]) ;
GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,0,[],...
    numChannelS,stdVal) ;
sigBinary = zeros(size(sigReshape)) ;
sigBinary(GammaBurstEvent.is_burst==1) = 1;
sigBinary = reshape(sigBinary,size(sigSmooth,1),size(sigSmooth,2),[]) ;

% prcPara = 95 ;
% prcBound = prctile(sigSmooth(:),prcPara) ;
% sigBinary = zeros(size(sigSmooth)) ;
% sigBinary(sigSmooth>prcBound) = 1;
%     
params.minPattTime = 30 ;
params.minPattSize = 3*resizeScale ;
flagMovie = 1 ;
flagSaveData = 1 ;
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;
saveFileName = [num2str(surMethodNum),'0.6SmRes8_3SD_',dataFileName,'_Band',...
   num2str(subBand(1)),'Hz',dateStr,'.mat'] ;
pattDetection_v4(sigSmooth, sigBinary, params, flagMovie, flagSaveData,saveFileName) ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagPostAnalysis
%% fig.2 in the paper is generated in this section
load('FullBurst2SD%ma027_032_Band30HzMay10_17:11.mat')
%% fig.2(A) snapshot of Gamma pattern
for iBurst = 1:length(patternIdx)
    [I,J,K] = ind2sub(size(sigSmooth),patternIdx{iBurst}) ;
    timeCount = 1;
    for iTime = min(K):max(K)
        binaryPattern{iBurst} = zeros(size(sigSmooth,1),size(sigSmooth,2),max(K)-min(K)+1) ;
        binaryPattern{iBurst}(I(K==iTime),J(K==iTime),timeCount) = 1 ;
        timeCount = timeCount+1 ;
    end
    instantPattern{iBurst} = binaryPattern{iBurst}.*sigSmooth(:,:,rangeFrame(iBurst,1):rangeFrame(iBurst,2)) ;
end

for iBurst = 21:41
    close all
    timeStep = 10 ;
    if size(instantPattern{iBurst},3)<timeStep*7+2
        continue
    end
    
    count = 1 ;
    figure;
    set(gcf,'Position',[680 401 1533 560]) ;
    for iTime = 2:timeStep:timeStep*7+2
        % iTime = rangeFrame(iBurst,1)+2:5:rangeFrame(iBurst,1)+5*7+2 % rangeFrame(iBurst,2)
        subplot(2,4,count)
        sigPlot = instantPattern{iBurst}(:,:,iTime) ;
        imagesc(sigPlot)
        set(gca,'YDir','normal')
%         title(['Detected burst at ',num2str(subBand),'Hz at ', ...
%             int2str(timeSlot),'ms'])
            caxis([0 0.5*max(sigSmooth(:))])
        colorbar
%         
%         imagesc(sigPlot(:,:,iTime))
%         set(gca,'YDir','normal')
%         % title(['95% signal at ',num2str(subBand),'Hz at ', ...
%         %    int2str(timeSlot/1000*1000),'ms'])
%         if plotAmp
%             caxis([0 0.5*max(sigSmooth(:))])
%         else
%             colorMapSpec = pmkmp_new;
%             sigLims = [-pi pi];
%             colormap(gca, colorMapSpec)
%             caxis(sigLims)
%         end
%         colorbar
        count = count+1 ;
        
    end
%    savefig(['PatternPropagation',num2str(iBurst),'.fig'])



%% fig.2(B) trajectory of Gamma pattern
    figure;
    set(gcf,'Position',[680 401 1533 560]) ;
    centre = WCentroids{iBurst} ;
    count = 1 ;
    for iTime = 2:1:5*7+2 % rangeFrame(iBurst,2)
        
        %subplot(2,4,count)
        plot([centre(iTime,1),centre(iTime+1,1)],[centre(iTime,2),...
            centre(iTime+1,2)],'k.-','markersize',20)

        count = count+1 ;
        xlim([0 10*resizeScale])
        ylim([0 10*resizeScale])
        hold on
    end
    pause
    close all
%    savefig(['PatternTrajectory',num2str(iBurst),'.fig'])
end
%% new pattern snapshots
load('FullBurst2SD_ma027_032_Band30HzJuly10_18:03.mat')
%%
sigTemp = zeros(size(sigSmooth)) ;
for iBurst = 1:length(patternIdx)
    sigTemp(patternIdx{iBurst}) = 1;
end

sigPattern = sigSmooth.*sigTemp ;
%%
for iBurst = 1:10
    centerCount = 1 ;
    for iTime = rangeFrame(iBurst,1):rangeFrame(iBurst,2)
        imagesc(sigPattern(:,:,iTime))
        hold on
        plot(WCentroids{iBurst}(centerCount,1),WCentroids{iBurst}(centerCount,2),'.', 'markersize',12)
        title(['burst number: ',num2str(iBurst)])
        pause
        cla
        centerCount = centerCount + 1 ;
    end
end

%% new pattern snapshots2'
burstNum = [1;2;3] ;
for iBur = 1:3
    iBurst = burstNum(iBur) ;
    subplot(3,2,(iBur-1)*2+1)
    iTime = rangeFrame(iBurst,1)+10 ;
    imagesc(sigPattern(:,:,iTime))
    hold on
    centre = WCentroids{iBurst} ;
    iTime = 10 ;
    plot([centre(iTime,1)],[centre(iTime,2),...
            ],'k.-','markersize',12)
%     xlim([0.5 20.5])
%     ylim([0.5 20.5])
    xlim([0.5 10.5])
    ylim([0.5 10.5])
    
    subplot(3,2,(iBur-1)*2+2)
    iTime = rangeFrame(iBurst,1)+40 ;
    imagesc(sigPattern(:,:,iTime))
    hold on
    
    centre = WCentroids{iBurst} ;
    count = 1 ;
    for iTime = 10:2:10+30 % rangeFrame(iBurst,2)
        if size(centre,1)<40
            continue
        end
        %subplot(2,4,count)
        plot([centre(iTime,1),centre(iTime+1,1)],[centre(iTime,2),...
            centre(iTime+1,2)],'k.-','markersize',8)

        count = count+1 ;
%         xlim([0.5 20.5])
%         ylim([0.5 20.5])
        xlim([0.5 10.5])
        ylim([0.5 10.5])

        hold on
    end
end

%% new pattern snapshots3'

lo = 0 ;
sdN2 = mean(sigSmooth(:))-2*std(sigSmooth(:)) ;
sdN1 = mean(sigSmooth(:))-1*std(sigSmooth(:)) ;
me = mean(sigSmooth(:));
sd1 = mean(sigSmooth(:))+1*std(sigSmooth(:)) ;
sd2 = mean(sigSmooth(:))+2*std(sigSmooth(:)) ;
sd3 = mean(sigSmooth(:))+3*std(sigSmooth(:)) ;
sd4 = mean(sigSmooth(:))+4*std(sigSmooth(:)) ;

%%
for burstNum = 1000:1600 ;    % 74,101,124; 171; 274
    % timePoint = [-5,5,15,25] ;
    timePoint = [-19,1,21,41] ;
    % if Duration(burstNum)<36
    if Duration(burstNum)<41+21+2
        continue
    end
    for iTime = 1:4
        timeFrame = rangeFrame(burstNum,1)+timePoint(iTime) ;
        iBurst = burstNum ;
        subplot(2,2,iTime)
        % imagesc(sigSmooth(:,:,timeRange))
        
        % hold on
        contourf(sigSmooth(:,:,timeFrame),[lo,sdN2,sdN1 ,me,sd1,sd2,sd3,sd4])
        
        if iTime >1
            hold on
            centre = WCentroids{iBurst} ;
            % iTimeS = (iTime-1)*10-5 ;
            iTimeS = (iTime-2)*20+1 ;
             plot([centre(iTimeS,1)],[centre(iTimeS,2)...
                    ],'k.','markersize',24,'linewidth',1)
            % iTimeS = (iTime-1)*10+4 ;
            iTimeE = (iTime-2)*20+21 ;
             plot([centre(iTimeE+2,1)],[centre(iTimeE+2,2)...
                    ],'ko','markersize',4,'linewidth',1)
                
            for iTime2 = iTimeS:2:iTimeE % rangeFrame(iBurst,2)
                
                %subplot(2,4,count)
                plot([centre(iTime2,1),centre(iTime2+2,1)],[centre(iTime2,2),...
                    centre(iTime2+2,2)],'r--','markersize',4,'linewidth',2)
                
                xlim([1 20])
                ylim([1 20])
                %         xlim([0.5 10.5])
                %         ylim([0.5 10.5])
                
                hold on
            end
        end
        
    end
    pause
    close all
end

%% 
thre2SD = mean(sigSmooth(:)) + 2*std(sigSmooth(:)) ;
thre1SD = mean(sigSmooth(:)) + 1*std(sigSmooth(:)) ;

for iBurst = 20:100
    centerCount = 1 ;
    for iTime = rangeFrame(iBurst,1):rangeFrame(iBurst,2)
        % imagesc(1:10,1:10,sigPattern(:,:,iTime))
        % hold on
        c = contourc(sigSmooth(:,:,iTime),[thre2SD,0]) ;
        plot((c(1,2:end)+0.5)/2,(c(2,2:end)+0.5)/2,'k-','markersize',12)
        xlim([0 10])
        ylim([0 10])
        %imagesc(sigPattern(:,:,iTime))
        hold on
        plot(WCentroids{iBurst}(centerCount,1)/2,WCentroids{iBurst}(centerCount,2)/2,'.', 'markersize',12)
        title(['burst number: ',num2str(iBurst)])
        pause
        cla
        centerCount = centerCount + 1 ;
    end
end

%% different colour of time
close all
colorStr = jet(8) ;
timeStep = 5 ;
timeStart = 3 ;
numStep = 8 ;
for iBurst = 627            %121 105 615 618 627 648 750  782 788 814 826 836 838 848 897 898 945 955 958 959 962
    centerCount = 1 ;
    for iTime = rangeFrame(iBurst,1)+timeStart+1:timeStep:...
            rangeFrame(iBurst,1)+timeStep*numStep+timeStart%rangeFrame(iBurst,2)
        % imagesc(1:10,1:10,sigPattern(:,:,iTime))
        % hold on
        
        c = contourc(sigSmooth(:,:,iTime),[thre2SD,0]) ;
        plot((c(1,2:end)+0.5)/2,(c(2,2:end)+0.5)/2,'k.--','MarkerEdgeColor',...
            colorStr(centerCount,:),'MarkerFaceColor',colorStr(centerCount,:),'markersize',12)
        xlim([0 10])
        ylim([0 10])
        %imagesc(sigPattern(:,:,iTime))
        hold on
        % plot(WCentroids{iBurst}(centerCount,1)/2,WCentroids{iBurst}(centerCount,2)/2,'.', 'markersize',12)
        
        centerCount = centerCount + 1 ;
    end
    % title(['burst number: ',num2str(iBurst)])
    legend('t = 0 ms','t = 5 ms','t = 15 ms','t = 20 ms','t = 25 ms','t = 30 ms','t = 35 ms','t = 40 ms')
    
    figure;
    colorStr2 = jet(40) ;
    centerCount = 1 ;
    for iTime2 = rangeFrame(iBurst,1)+timeStart+1: rangeFrame(iBurst,1)+timeStep*8+timeStart
        plot(WCentroids{iBurst}(iTime2-rangeFrame(iBurst,1),1)/2,...
            WCentroids{iBurst}(iTime2-rangeFrame(iBurst,1),2)/2,...
            '.-', 'markersize',12,'MarkerEdgeColor',colorStr2(centerCount,:))
        xlim([0 10])
        ylim([0 10])
        centerCount = centerCount + 1 ;
        hold on
    end
    colorbar; caxis([0 40]);colormap(jet(40))
    %pause
    %cla
end






%% fig.2(C) MSD and distribution
center = WCentroids ;
% % select bursts of duration more than 100ms
% center = [] ;
% count = 1 ;
% for iBurst = 1:length(Duration)
%     if Duration(iBurst)>100
%         center{count} = WCentroids{iBurst} ;
%         count = count+1 ;
%     end
% end

fullMSD = [] ;
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 1 ;
slopeAll = [] ;
superIdx = [] ;
% close all
for iBurst = 1:size(center,2)     % 198
    Trajectory = [] ;
   Trajectory = center{iBurst} ;
   Trajectory(:,3) = 1:size(center{iBurst},1) ;  
    
% for iCentre = 1:floor(size(center{iBurst} ,1)/4)
%     Trajectory(iCentre,:) = mean(center{iBurst}((iCentre-1)*4+1:4*iCentre,:)) ;
% end
% if length(1:4:size(center{iBurst},1)) > size(Trajectory,1)
%     Trajectory(:,3) = 1:4:size(center{iBurst},1)-4 ;
% else
%     Trajectory(:,3) = 1:4:size(center{iBurst},1) ;
% end
    itaQ = nan(5,1) ;
    countQ = 1 ;
    
    for powerQ = 2% 1:0.2:5
        % [MSD,tau] = get_MD(Trajectory,powerQ) ;%,grid_size);
        [MSD,tau] = get_MSD(Trajectory) ;%,grid_size);
    tempMSD = [MSD,tau] ;
    fullMSD = [fullMSD;tempMSD] ;
    
    p_temp = [];
    norErr = [] ;
    maxStart = 30 ;               % minimum tau used to fit MSD
    tau_max = maxStart:min(100,size(MSD,1)) ;
    
    for fitIdx = 1:length(tau_max)       
        fitRange = (1:tau_max(fitIdx)) ;
        [pAll,S] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
        y = exp(polyval(pAll,log(tau(fitRange)))) ;
        errorRate = mean(abs((MSD(fitRange)-y(fitRange))./MSD(fitRange))) ;
        p_temp(fitIdx,:) = pAll ;
        % norErr(fitIdx) = S.normr ;       
        norErr(fitIdx) = errorRate ;
    end
    % discard the MSD with error rate more than 10%
    if (min(norErr)<0.1)
        bestIdx = find(norErr == min(norErr)) ;
        pBest = p_temp(bestIdx(end),:) ;
        y = exp(polyval(pBest,log(tau))) ;
        slopeAll(iBurst) = pBest(1) ;
    else
       continue
    end
    % keep the superdiffusive one to do average    
    %if pBest(1)>1
        msdSuper = [msdSuper;tempMSD] ;
        superIdx(countSuper) = iBurst ;
        countSuper = countSuper + 1; 
    %end
    %plot for single burst
%     loglog(tau,MSD,'.-','MarkerSize',12)
%     % plot(tau,MSD,'.-')
%     
%     hold on
%     loglog(tau,y)
%     % plot(tau,y)
%     title(['MSD within a burst (',num2str(size(center{iBurst},1)),' ms) '])%,...
%         %'with \alpha = ', num2str(pAll(1)),' tau = ',num2str(maxStart+bestIdx(end))])
%     xlabel('\tau (ms)')
%     ylabel('Mean Square Distance (mm)')
%     
%     str = {'p = ',num2str(pAll(1))};
%     text(max(tau)+5,max(y)+5,str)
%     xlim([1 10^3])
%     pause
%     close all
    % itaQ(countQ) = pAll(1) ;
    % countQ = countQ+1 ;
    end
    % plot(linspace(1,5,length(itaQ)),itaQ,'ro-')
    % pause
    % close all
end

% plot for average MSD and distribution
figure;
sigIn = msdSuper ;
[MSDsort,~,MSDIdx] = unique(sigIn(:,2)) ;

for idx = 1:length(MSDsort)
    meanMSD(idx) = mean(sigIn(find(MSDIdx==idx),1)) ;
    stdMSD(idx) = std(sigIn(find(MSDIdx==idx),1)) ;
end

% plot(MSDsort,meanMSD,'o')
% errorbar(MSDsort,meanMSD,stdMSD)
loglog(MSDsort,meanMSD)
% set(gca,'yscale','log','xscale','log')
hold on

tau = 1:200 ;
p = polyfit(log(MSDsort(1:20)'),log(meanMSD(1:20)),1) ;
y = exp(polyval(p,log(tau))) ;
slope = p(1) ;
loglog(tau,y)
% title(['Mean Square Distance versus \tau within a burst (',num2str(size(center{iBurst},1)),')'])
xlabel('\tau (ms)')
ylabel('Mean Square Distance (electrode)')

str = {'p = ',num2str(p(1))};
text(max(tau)+5,max(y)+5,str)

figure
slopeAll(slopeAll == 0) = [] ;
hist(slopeAll,40)
% title('distribution of MSD for my144')
xlabel('MSD')


% fig.2(D) collapsed displacement distribution
%% displacement distribution
figure;
center = WCentroids ;
deltaTime = [2,4,8,16,32] ; %,100]; %30 ;   
iTa = 0.64 ;
% iTa = 0.64 ;

legendInfo = [] ;
count = 0 ;
for delta = 1:length(deltaTime) %:-1:1
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    DisplaceNorm = [] ;
for iBurst = 1:length(superIdx) %  size(center,2)  % 
    curBurst = superIdx(iBurst) ;
    % curBurst = iBurst ;
    posCenter = center{curBurst} ;
    % DisplaceTemp = sqrt(sum((posCenter(deltaT+1 :end,:) - posCenter(1: end-deltaT,:)).^2,2)) ;
    DisplaceTemp = (posCenter(deltaT+1 :end,1) - posCenter(1: end-deltaT,1) ) ;
    Displace = [Displace;DisplaceTemp] ;
    DisplaceTempNorm = DisplaceTemp/(deltaT^iTa) ;
    DisplaceNorm = [DisplaceNorm;DisplaceTempNorm] ;
end
count = count+1 ;
if count>1
    allDisplaceTemp = [DisplaceNorm; allDisplaceTemp] ;
else
    allDisplaceTemp = DisplaceNorm;
end

[n,x] = histcounts(DisplaceNorm,'normalization','pdf') ;
% loglog(x(2:end),n,'o')
semilogy(0.5*(x(1:end-1)+x(2:end)),n,'o')
hold on
legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
end

fitVar = allDisplaceTemp ; % allDisplaceTemp ;  
pd = fitdist(fitVar,'stable') ;
x = linspace(x(1)-3,x(end)+3,100000) ;
y = pdf(pd,x) ;
hold on ; % loglog(x,y,'r-');
semilogy(x,y,'r-');

pd = fitdist(fitVar,'normal') ;
y = pdf(pd,x) ;
hold on ; % loglog(x,y,'k--');
semilogy(x,y,'k--');

legend('t = 30 ms','Gaussian')
legendInfo{delta+1} = ['\alpha stable fit'] ;
legendInfo{delta+2} = ['Gaussian fit'] ;
legend(legendInfo)
xlabel('electrode')
% ylabel('displacement/\eta^{0.64}')
% title('Displacement distribution')


%%
makeMovie = 1;
if makeMovie
close all
timeUnit = 's';
t = datetime('now') ;
dateStr = datestr(t,'mmmmdd_HH:MM') ;
vidTitle = [pwd,'/Results/Project1/',dataFileName,'_',...
        dateStr,'.avi'] ;
    vidObj = VideoWriter(vidTitle,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 20 ;
    open(vidObj);
    fig=figure ;
    % set(gcf,'Visible','On','Position',[104 443 1880 455])
    set(gcf,'Visible','On','Position',[31 299 1525 587])
    set(gcf,'Visible','On','Position',[360 30 790 665])
    
     
for iTime = fix(30*fsTemporal)+1:1:fix(40*fsTemporal)  %fix(1.08*fsDS):fix(1.16*fsDS)
        subplot(2,2,1)
    % % Set signal and colors based on whether you are plotting amplitude
    % Amp
    % Optionally interpolate spatial grid of data
    signalGrid = abs( squeeze(sigIn(:,:,iTime)) ) ;
    
    colorMapSpec = parula;
    sigLims = [min(abs(sigIn(:))), 0.4*max(abs(sigIn(:)))];
    % sigLims = [1e-6 8e-6] ;
    imagesc(linspace(1,10,40),linspace(1,10,40),signalGrid,sigLims)
    title(['original Gamma at ',num2str(iTime/fsTemporal,'%2.3f'), 's' ])
    
    
    subplot(2,2,2)
    signalGrid = abs( squeeze(sigSmooth(:,:,iTime)) ) ;
    colorMapSpec = parula;
    sigLims = [min(abs(sigSmooth(:))), 0.4*max(abs(sigSmooth(:)))];
    imagesc(linspace(1,10,40),linspace(1,10,40),signalGrid,sigLims)
     title(['Smoothed Gamma amplitude at ',num2str(iTime/fsTemporal,'%2.3f'), 's' ])
 
         subplot(2,2,3)
    signalGrid = abs( squeeze(sigPlot(:,:,iTime)) ) ;
    colorMapSpec = parula;
    sigLims = [min(abs(sigPlot(:))), 0.4*max(sigPlot(:))];
        imagesc(linspace(1,10,40),linspace(1,10,40),signalGrid,sigLims)
     title(['Gamma burst detection at ',num2str(iTime/fsTemporal,'%2.3f'), 's' ])
     
         subplot(2,2,4)
    signalGrid = abs( squeeze(patternPlot(:,:,iTime)) ) ;
    colorMapSpec = parula;
    sigLims = [min(abs(patternPlot(:))), 0.4*max(abs(patternPlot(:)))];
        imagesc(linspace(1,10,40),linspace(1,10,40),signalGrid,sigLims)
     title(['Gamma patterns at ',num2str(iTime/fsTemporal,'%2.3f'), 's' ])
     


    writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    cla
end
close(vidObj);

end


%% suplementary for finding peak frequency distribution
close all
tempData = squeeze (sigOri(5,5,:)) ;
[wt2,f22,coi] = cwt(tempData,fsTemporal,'VoicesPerOctave',30);

temp = abs(wt2(75:118,:));
temp2 = prctile(temp,95) ;
temp3 = zscore(temp,[],2) ;
temp3 = temp ;
temp4 = prctile(temp3(:),95) ;
temp5 = zeros(size(temp)) ;
temp5(temp3>temp4) = 1;

CC = bwconncomp(temp5) ;    
for iBurst = 1:length(CC.PixelIdxList)
    [~,idx] =max(temp3(CC.PixelIdxList{iBurst})) ;
    peakFreqTemp = mod(CC.PixelIdxList{iBurst}(idx),size(temp5,1)) ;
    if peakFreqTemp == 0
        peakFreqTemp = size(temp5,1) ;
    end
    peakFreq(iBurst) = f22(peakFreqTemp+75-1) ;
end

histogram(peakFreq,30)


%% suplementary for alpha stable fit
sigDist = cell2mat(GammaBurstEvent.burst_du_steps) ;
% sigDist = cell2mat(GammaBurstEvent.flat_du_steps) ;
paraDist = fitdist(sigDist','Stable') ;
[x,n] = histcounts(sigDist,'normalization','pdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = pdf(paraDist,n0) ;
figure;
semilogy(n0,x,'o')
hold on
semilogy(n0,y0,'r')

%%
sigDist = cell2mat(GammaBurstEvent.burst_du_steps) ;
% sigDist = cell2mat(GammaBurstEvent.flat_du_steps) ;
paraDist = fitdist(sigDist','LogNormal') ;
[x,n] = histcounts(sigDist,'normalization','pdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = pdf(paraDist,n0) ;
figure;
semilogy(n0,x,'o')
hold on
semilogy(n0,y0,'r')

%%
load('FullBurst2SD_ma027_032_Band30HzJuly10_18:03.mat')
sigDist = pattSize ;
sigDist = Duration ;
sigDist = centInterval ;

%%
sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
%paraDist = fitdist(sigDist','Stable') 
[x,n] = histcounts(sigDist,200,'normalization','pdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = 2*pdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
figure;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
semilogy(n0,x,'b.','MarkerSize',12)
hold on
semilogy(n0,y0,'r')
legend('original data','\alpha stable fit')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

[R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable(sigDist) ;

%% binplot
sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
%paraDist = fitdist(sigDist','Stable') 
figure;
numBins = 45 ;
histogram(sigDist,numBins,'normalization','pdf') ;
% set(gca,'XScale','log')
set(gca,'YScale','log')
[x,n] = histcounts(sigDist,numBins,'normalization','pdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = 2*pdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
% semilogy(n0,x,'b.','MarkerSize',12)
hold on
 semilogy(n0,y0,'r')
% loglog(n0,y0,'r')

legend('original data','\alpha stable fit')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

[R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable(sigDist) ;

%% cdf
sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
% paraDist = fitdist(sigDist','Stable') 
%paraDist = fitdist(sigDist','Stable') 
[x,n] = histcounts(sigDist,100,'normalization','cdf') ;
n0 = 0.5*(n(2:end)+n(1:end-1)) ;
y0 = cdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
figure;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
loglog(n0,1-x,'b.','MarkerSize',12)
hold on
loglog(n0,1-y0,'r')
legend('original data','\alpha stable fit')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

%%
sigDist(sigDist<=0) = [] ;
sigFull = [sigDist,-sigDist] ;
paraDist = fitdist(sigFull','Stable') 
%paraDist = fitdist(sigDist','Stable') 
nEdge = logspace(log10(min(sigDist)),log10(max(sigDist)),20) ;

[x,n] = histcounts(sigDist,nEdge,'normalization','pdf') ;
n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ;
y0 = 2*pdf(paraDist,n0) ;
%y0 = pdf(paraDist,n0) ;
figure;
% loglog(n0,x,'b.','MarkerSize',12)
% hold on
% loglog(n0,y0,'r')
loglog(n0,x,'b.','MarkerSize',12)
hold on
loglog(n0,y0,'r')
legend('original data','\alpha stable fit')
% xlabel('size (sites)')
% xlabel('duration (ms)')
xlabel('intervel (ms)')
ylabel('probability')

[R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable(sigDist) ;
end

