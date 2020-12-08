function hilbertSig = find_Hilbert(sigIn, fsTemporal, timeDim, badChannels,...
    discardTime)

% hilbertSig = find_Hilbert(sigIn, fsTemporal,...
%    badChannels, discardTime)
%
% This function finds the bandpass signal using a butterworth two way
% filtering
%
% Input: 
% sigIn            (1 X Time) or (Freq X ChannelX X ChannelY X Time) signal
% fsTemporal       sampling frequency in time
% badChannels      dead electrodes
% discardTime      discard time (in second) for filtering, default = 
%                  min (0.2s, 1/50 * total timelength)
%
% Output:
% hilbertSig      (Freq X ChannelX X ChannelY X Time) bandpass signals

% author: Xian Long   Supervisor: Pulin Gong
% Nov, 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default settings
if ~exist('badChannels', 'var') 
    switch timeDim
        case 4
            nanChans = any(isnan(squeeze(sigIn(1,:,:,:))),timeDim-1) ;
            zeroChans = all(squeeze(sigIn(1,:,:,:))==0, timeDim-1) ;
            badChannels = find(nanChans | zeroChans) ;
    end
end
if ~exist('discardTime', 'var')
    discardTime = min( 0.2, ( 0.02*size(sigIn,timeDim)/fsTemporal ) ) ;
end

if fix(size(sigIn,timeDim)/fsTemporal) < discardTime*2
    error('The length of signals is smaller than discard time')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
if timeDim == 4
    gridLength = size(sigIn, 2) ;
elseif timeDim == 2
    gridLength = 1 ;
else
    error('Raw signal format incorrect!') 
end

reshapeSig = reshape(sigIn, size(sigIn,1), gridLength*gridLength,[]) ;
hilbertFull = nan(size(reshapeSig));
discardTimeSteps = fix(discardTime * fsTemporal);
outSection = discardTimeSteps+1 : (size(reshapeSig,3)-discardTimeSteps) ;

for bandNum = 1:size(sigIn,1)
    % Filter LFPs
    disp('Hilbert Transform...') ; 
       
    for ichannel = setdiff(1:size(reshapeSig,2), badChannels)    
        hilbertFull(bandNum,ichannel,:) = hilbert(reshapeSig(bandNum,ichannel,:));
    end
    disp(['percentage ',int2str(bandNum/size(sigIn,1)*100),'%'])
end
% Cut off the start and ends of the signal that may experience ...
% boundary effects from the filtering signals
hilbertSig = hilbertFull(:,:,outSection) ;

hilbertSig = reshape(hilbertSig,size(hilbertSig,1),gridLength,...
    gridLength,length(outSection)) ;
