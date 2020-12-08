function postAnalysisP2(sigIn,saveFolderName,fsTemporal)

%% correlation decay between frequencies
disp(['finding correlation decay'])
bandName = [{'1Hz'},{'4Hz'},{'16Hz'},{'64Hz'}] ;
legendName = [] ;
% close all
 figure;
 Band = [1,5,9,13] ;
% Band = [5,6,7] ;
for iBandIdx = [1:4]
    iBand = Band(iBandIdx) ;
timeLength = fix(30*fsTemporal) ;
timeStart = 1 ;
sigCor = squeeze(sigIn{iBand}(:,:,:)) ;

    timeSeries = sigCor(:,:,timeStart:timeStart+timeLength) ;
    timeSeries = reshape(timeSeries,size(timeSeries,1)*size(timeSeries,2),[]) ;
    timeSeries = timeSeries' ;
    R = corr(timeSeries) ;

% find the distance
d = [] ;
D = [] ;
 [X,Y] = meshgrid(1:size(sigCor,1),1:size(sigCor,2)) ;
x = X*0.2 ;   % 0.4mm electrode spacing
y = Y*0.2 ;   % 0.4mm electrode spacing

for ix = 1:size(X,1)*size(X,2)
    X1 = [X(ix),Y(ix)] ;
    x1 = [x(ix),y(ix)] ;
    for iy = 1:size(X,2)*size(X,1)
        X2 = [X(iy),Y(iy)] ;
        x2 = [x(iy),y(iy)] ;
        d(ix,iy) = sqrt(sum((x1-x2).^2)) ;
        D(ix,iy) = sqrt(sum((X1-X2).^2)) ;
    end
end
% find the average correlation
Rsort = [] ;
 [C,IA,IC] =unique(d) ;
 for iDist = 1:length(C)
     idx = find(d==C(iDist)) ;
     Rsort(iDist) = nanmean(R(idx)) ; 
 end

 plot(C,Rsort','o')   %.*(subBand(iBand).^0.4)
legendName{iBandIdx} = [bandName{iBandIdx}] ;
hold on
end
legend(legendName)
title('correlation across space over different frequency bands (50s signals)')
grid on
xlabel('distance (mm)')
ylabel('correlation coef.')

% saveas(gcf,[pwd,'/Results/Project2/final/correlationLength/',dataFileName,dateStr,'.eps'],'epsc') ;

%% scaling behaviour of pattern properties in different frequencies
disp(['finding scaling between freq'])
% Find subband signals
% apply wavelet filtering for Gamma signals
% numBand = 15 ;
% wvcfs = cell(numBand,1) ;
% cfreqso = zeros(numBand,1) ;
% 
% for iBand = 1:numBand
%     cfreqso(iBand) = 2^((iBand-1)/2) ;
%     [wvcfs{iBand},~] = find_cwtCoef(sigSmooth,cfreqso(iBand), fsTemporal,7) ;
% end
% 
% for iBand = 1:numBand
%     sigIn{iBand} = abs(wvcfs{iBand}) ;
% end

%% scaling plot

scaleBand = zeros(numBand,1) ;
instScaleBand = zeros(numBand,1) ;
for iNum = 1:numBand
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
    load([saveFolderName,fileName(1).name],'patternScale')
    load([saveFolderName,fileName(1).name],'instantScale')
    scaleBand(iNum) = mean(patternScale) ;
    for iBurst = 2:length(instantScale)
        maxSize(iBurst) = max(instantScale{iBurst}) ;
    end
    instScaleBand(iNum) = mean(maxSize) ;
end

for iNum = 1:numBand
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','my144*']) ;
    load([saveFolderName,fileName(1).name],'patternScale')
    load([saveFolderName,fileName(1).name],'instantScale')
    scaleBand(iNum) = mean(patternScale) ;
    for iBurst = 2:length(instantScale)
        maxSize(iBurst) = max(instantScale{iBurst}) ;
    end
    instScaleBand(iNum) = mean(maxSize) +instScaleBand(iNum) ;
end

for iNum = 1:numBand
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','my147*']) ;
    load([saveFolderName,fileName(1).name],'patternScale')
    load([saveFolderName,fileName(1).name],'instantScale')
    scaleBand(iNum) = mean(patternScale) ;
    for iBurst = 2:length(instantScale)
        maxSize(iBurst) = max(instantScale{iBurst}) ;
    end
    instScaleBand(iNum) = mean(maxSize) +instScaleBand(iNum) ;
end

% for iNum = 1:numBand
%     maxSize = [] ;
%     
%     fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma025*']) ;
%     load([saveFolderName,fileName(1).name],'patternScale')
%     load([saveFolderName,fileName(1).name],'instantScale')
%     scaleBand(iNum) = mean(patternScale) ;
%     for iBurst = 2:length(instantScale)
%         maxSize(iBurst) = max(instantScale{iBurst}) ;
%     end
%     instScaleBand(iNum) = mean(maxSize) +instScaleBand(iNum) ;
% end

for iNum = 1:numBand
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
    load([saveFolderName,fileName(2).name],'patternScale')
    load([saveFolderName,fileName(2).name],'instantScale')
    scaleBandSur(iNum) = mean(patternScale) ;
    for iBurst = 2:length(instantScale)
        maxSizeSur(iBurst) = max(instantScale{iBurst}) ;
    end
    instScaleBandSur(iNum) = mean(maxSizeSur) ;
end

for iNum = 1:numBand
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma025*']) ;
    load([saveFolderName,fileName(2).name],'patternScale')
    load([saveFolderName,fileName(2).name],'instantScale')
    scaleBandSur(iNum) = mean(patternScale) ;
    for iBurst = 2:length(instantScale)
        maxSizeSur(iBurst) = max(instantScale{iBurst}) ;
    end
    instScaleBandSur(iNum) = mean(maxSizeSur)+ instScaleBandSur(iNum);
end


figure;
loglog(cfreqso,scaleBand,'o')
hold on
loglog(cfreqso,scaleBandSur,'r--')

figure;
loglog(cfreqso,instScaleBand/3,'o')
hold on
loglog(cfreqso,instScaleBandSur/2,'r--')
ylim([10,120])

%% scaling behaviour of pattern properties within frequencies
disp(['finding scaling within freq'])
figure;
% saveFolderName = [pwd,'/Results_data/Project2/BurstDetection3/'] ;
addpath(genpath([pwd,'/ToolOthers/Clauset_powerlaws/']))
scaleBand = zeros(15,1) ;
instScaleBand = zeros(15,1) ;

iNum = 9 ;
fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
load([saveFolderName,fileName(1).name],'Duration')
load([saveFolderName,fileName(1).name],'patternScale')
load([saveFolderName,fileName(1).name],'instantScale')

sigIn = patternScale ;
% sigIn = Duration ;

sigMinScale = mean(sigIn) ;
cfreqMinScale = cfreqso(iNum) ;

for iNum = 9:2:17
    maxSize = [] ;
    
    fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
    load([saveFolderName,fileName(1).name],'Duration')
    load([saveFolderName,fileName(1).name],'patternScale')
    load([saveFolderName,fileName(1).name],'instantScale')
    
    
    sigIn = patternScale ;
%     fileName = dir([saveFolderName,num2str(iNum,'%02d'),'FullBurst95%*','ma027*']) ;
%     load([saveFolderName,fileName(2).name],'Duration')
%     load([saveFolderName,fileName(2).name],'patternScale')
%     load([saveFolderName,fileName(2).name],'instantScale')
%   sigIn = Duration ;
    
    xmin =  min(sigIn(:)) ;
    % sigIn = sigIn(sigIn>xmin) ;
    
    pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
    [lambdaHat1,lambdaCI] = mle(sigIn, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
    L1 = sum(log(pf_power(sigIn,lambdaHat1 ))) ;
    
    numPoints = 200 ;
    z = linspace(min(sigIn(:)),max(sigIn(:)),numPoints);
    x = 0.5*(z(1:end-1)+z(2:end));
    % histogram(sigIn,z,'Normalization','pdf');
    
    % hold on
    % plot(x,pf_power(x,lambdaHat1))
    
%     figure
    % nEdge = logspace(log10(5000),log10(max(sigIn(:))),20) ;

    [c,n] = histcounts(sigIn,z,'Normalization','pdf');
    % n0 = 10.^(0.5*(log10(n(2:end))+log10(n(1:end-1)))) ;

    gamma = 0.24^(cfreqso(iNum)/cfreqso(17)) ;
    gamma = 1/cfreqso(iNum)^2*sigMinScale* cfreqMinScale^2 /mean(sigIn) ;
    gamma = 1 ;

    loglog(x/gamma ,c,'o','markersize',4,'linewidth',2)
     hold on
%     loglog(x,pf_power(x, lambdaHat1),'k');
    % surrogate data
    % hold on
%     z = linspace(min(sigSur(:)),max(sigSur(:)),numPoints);
%     x = 0.5*(z(1:end-1)+z(2:end));
%     [c2] = histcounts(sigSur,z,'Normalization','pdf');
%     loglog(x,c2,'o','markersize',4,'linewidth',2)
%     hold on
    % MLE of exponential
% hold on;
% pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
% [lambdaHat3,lambdaCI] = mle(sigSur, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
% loglog(x,pf_exp(x, lambdaHat3),'r');
    
%     % MLE of power law with cut off
%     figure;
%     pf_powerC = @(x,alpha,lambda) x.^-alpha.*exp(-lambda*x)*lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
%     [lambdaHat2,lambdaCI] = mle(sigIn, 'pdf',pf_powerC, 'start',[1,0], 'lowerbound',[1,0], 'upperbound', [3,10]) ;
%     L2 = sum(log(pf_powerC(sigIn,lambdaHat2(1), lambdaHat2(2) ))) ;
%     
%     histogram(sigIn,z,'Normalization','pdf');
%     
%     hold on
%     plot(x,pf_powerC(x,lambdaHat2(1),lambdaHat2(2)))
    
    
end

% powerLawSup