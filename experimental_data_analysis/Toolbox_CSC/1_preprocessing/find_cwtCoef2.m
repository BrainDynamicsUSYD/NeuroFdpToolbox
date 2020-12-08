function [cwtCoef,pseudoFreq] = find_cwtCoef2(sigIn, freqRange, scaleStep,...
    fsTemporal,wName, badChannels)

% [cwtCoef,pseudoFreq] = find_cwtCoef2(sigIn, freqRange, scaleStep,...
%    fsTemporal,wName) ;
%
% This function finds the complex wavelet coefficient by calling the cwt
% function.
% 
% Inputs:
%       sigIn          (1 X Time) or (ChannelX X ChannelY X Time) raw
%                           signals
%       freqRange      [lowFreqLimit highFreqLimit]
%       scaleStep      scale step in cwt, default = 1;
%       fsTemporal     sample frequency
%       wname          wavelet name, default = 'cmorl1.5-1'
%       badChannels    bad channels
%
% Outputs:
%       cwtCoef        (Freq X Time) or (Freq X ChannelX X ChannelY X Time)
%                          Complex Wavelet Coefficient
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings
if ~exist('badChannels', 'var') 
    nanChans = any(isnan(sigIn(:,:,:)),3);
    zeroChans = all(sigIn(:,:,:)==0, 3);
    badChannels = find(nanChans | zeroChans);
end
if ~exist('scaleStep','var')
    scaleStep = 1;
end
if ~exist('wName','var')
    wName = 'cmor1.5-1';             % cmor fc-fb tune this to have a proper
                                     % resolution
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
gridLength = size(sigIn, 1) ;
reshapeSig = reshape(sigIn, gridLength*gridLength,[]) ;

% wavelet scales   (standard steps from the matlab example)
fc = centfrq(wName) ;
scalerange = fc./(freqRange/fsTemporal) ;
scales = scalerange(end):scaleStep:scalerange(1) ;            % 
pseudoFreq = scal2frq(scales, wName, 1/fsTemporal) ;

% loop through each channel to find wavelet coeficient
cwtCoefTemp = zeros(length(pseudoFreq),gridLength*gridLength,size(reshapeSig,2)) ;
for ichannel = setdiff(1:size(reshapeSig,1), badChannels)
    cwtCoefTemp(:,ichannel,:) = cwt( reshapeSig(ichannel,:) ,scales, wName  ) ;
    disp(['percentage ',int2str(ichannel/size(reshapeSig,1)*100),'%'])
end

% reshape the data to the input format
cwtCoef = squeeze(reshape(cwtCoefTemp,size(cwtCoefTemp,1),gridLength,...
    gridLength,[]))  ;
