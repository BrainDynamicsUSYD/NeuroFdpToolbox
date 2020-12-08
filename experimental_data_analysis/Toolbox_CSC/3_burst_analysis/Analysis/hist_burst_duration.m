%% Gamma burst duration distribution without any threshold (95%)
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
subBand = [30,80] ;
bandpassSig2 = find_bandpassSig(sigOri,subBand, fsTemporal,3,badChannels,0) ;
%%
hilbertSig2 = find_Hilbert(bandpassSig2, fsTemporal,4,badChannels) ;
sigIn2 = abs(squeeze(hilbertSig2(1,:,:,fix(11*fsTemporal)+1:fix(100*fsTemporal)) ));

%%
boundPrctile = 95 ;
prcBound = prctile(sigIn2(:),boundPrctile) ;
sigBinary = zeros(size(sigIn2)) ;
sigBinary(sigIn2>prcBound) = 1;
sigPlot = sigBinary.*sigIn2 ;

%%
distDuration = [] ;
for xPoint = 1:10
    for yPoint = 1:10
        tempSig = sigBinary(xPoint,yPoint,:) ;
        tempS = tempSig(1:end-1) - tempSig(2:end) ;
        burstStart = find(tempS == 1) ;
        burstEnd = find(tempS == -1) ;
        if isempty(burstStart)||isempty(burstEnd)
            continue
        end
        if burstStart(1)<burstEnd(1)
            if  burstStart(end)<burstEnd(end)
                duration = burstEnd - burstStart;
            else
                duration = burstEnd - burstStart(1:end-1) +1;
                duration = [duration; length(tempS)- burstStart(end)+1] ;
            end
        else
            if  burstStart(end)<burstEnd(end)
                duration = burstEnd(2:end) - burstStart + 1;
                duration = [burstEnd(1) ; duration] ;
            else
                duration = burstEnd(2:end) - burstStart(1:end-1);
                duration = [burstEnd(1) ; duration; length(tempS)- burstStart(end)+1] ;
            end
        end
        allDuration{xPoint,yPoint} = duration ;
        distDuration = [distDuration;duration] ;
    end
end

