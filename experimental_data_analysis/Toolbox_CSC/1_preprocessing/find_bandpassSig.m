function bandpassSig = find_bandpassSig(sigIn,subBand, fsTemporal,timeDim,...
    badChannels, normalFlag, filterOrder, discardTime)

%bandpassSig = find_bandpassSig(sigIn,subBand, fsTemporal,...
%    badChannels,gridLength, normalFlag, filterOrder, discardTime) ;
%
% This function finds the bandpass signal using a butterworth two way
% filtering
%
% Inputs: 
%       sigIn            (1 X Time) or (ChannelX X ChannelY X Time) raw
%                           signals
%       subBand          (FreqBandNums X 2) frequency range of the pass band
%       fsTemporal       sampling frequency in time
%       timeDim          the dimension of time, 2 or 3.
%       badChannels      dead electrodes
%       normalFlag       boolean variable, default = 0 for NO zscore
%       filterOrder      default = 8
%       discardTime      discard time (in second) for filtering, default = 
%                           min(1s, 1/10 * total timelength)                  
%
% Outputs:
%       bandpassSig      (FreqBandNum X Time) or (FreqBandNum X ChannelX X 
%                           ChannelY X Time) bandpass signals 
%
% author: Xian Long   Supervisor: Pulin Gong
% Nov, 2017

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default settings
if ~exist('badChannels', 'var') 
    nanChans = any(isnan(sigIn(:,:,:)),3);
    zeroChans = all(sigIn(:,:,:)==0, 3);
    badChannels = find(nanChans | zeroChans);
end
if ~exist('normalFlag', 'var') 
    normalFlag = 0 ;
end
if ~exist('filterOrder', 'var')
    filterOrder = 8 ;
end
if ~exist('discardTime', 'var')
    discardTime = min( 1, ( 0.1*size(sigIn,timeDim)/fsTemporal ) ) ;
end

if fix(size(sigIn)/fsTemporal) < discardTime*2
    error('The length of signals is smaller than discard time')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
gridLength = size(sigIn, 1) ;
reshapeSig = reshape(sigIn, gridLength*gridLength,[]) ;
discardTimeSteps = fix(discardTime * fsTemporal);
outSection = discardTimeSteps+1 : (size(reshapeSig,2)-discardTimeSteps);
bandpassSig = nan(size(subBand,1), gridLength*gridLength, length(outSection) ) ;

% perform bandpass filtering for each band defined
for bandNum = 1:size(subBand,1)
    % Filter Signals
    disp('Filtering waveforms...') ; 
    fLow = subBand(bandNum,1) ;
    fHigh = subBand(bandNum,2) ;
    % Design Butterworth band-pass filter
    h = fdesign.bandpass('N,F3dB1,F3dB2',filterOrder,fLow,fHigh,fsTemporal);
    Hd = design(h, 'butter');
    set(Hd, 'Arithmetic', 'double');
    
    SOS = Hd.sosMatrix;
    G = Hd.ScaleValues;
    bandpassSigFull = nan(size(reshapeSig));
    for ichannel = setdiff(1:size(reshapeSig,1), badChannels)
        % Filter signals forwards and backwards to avoid phase distortion
        ifx = filtfilt(SOS,G,reshapeSig(ichannel,:));
        bandpassSigFull(ichannel,:) = ifx;
        if normalFlag
            % normalise bands
            bandpassSigFull(ichannel,:) = zscore(bandpassSigFull(ichannel,:)) ;
        end
    end

    % Cut off the start and ends of the signal that may experience ...
    % boundary effects from the filtering signals    
    bandpassSig(bandNum,:,:) = bandpassSigFull(:,outSection) ;
    disp(['percentage ',int2str(bandNum/size(subBand,1)*100),'%'])
end

% reshape data to the input format
bandpassSig = (reshape(bandpassSig,size(bandpassSig,1),gridLength,...
    gridLength,length(outSection)))  ;