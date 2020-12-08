function [dataOut,fsTemporal, badChannels] = preprocess_LFP(dataIn, flagBandstop)

% [dataOut,fsTemporal,badChannels] = preprocess_LFP(dataIn, flagBandstop) ;
%
% This function implements the procdure to load the local field potential 
% data (spatial re-arrangement) from Prof. Paul Martin's group and includes
% a bandstop filter at [49,51Hz]. 
% 
% Inputs:
%       dataIn         a (Channel X Time) matrix of raw LFP data
%       flagBandstop   flag for bandstoping electric line, default = 1;
%
% Outputs:
%       dataOut        a (ChannelX X ChannelY X Time) matrix of LFP data
%       fsTemporal     sampling frequnecy
%       badChannels    channels that are NaN or zeros.
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings
if ~exist('flagBandstop','var')
    flagBandstop = 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Sampling rate in time
fsTemporal = 1.017252624511719e+03 ;
% Convert to X-Y 10*10*T*Trial grid of electrodes
dataOut = vector2grid(dataIn) ;

% Find any channels with NaN or zero values
nanChans = any(isnan(dataOut(:,:,:)),3) ;
zeroChans = sum(dataOut(:,:,:),3) == 0 ;
badChannels = find(nanChans | zeroChans);
% Also exclude corner electrodes from Utah MEA LFPs
badChannels = union([1 10 91 100], badChannels) ;

% bandstop the electric line (49-51Hz), use 6th order butterworth filter
if flagBandstop
    % initialisation
    reshapeLFPs = reshape(dataOut, 10*10,[]) ;
    filteredLFPs = nan(size(reshapeLFPs));
    
    % bandstop the 49-51 Hz, use 2nd order butterworth filter
    fLowNorm = 49/fsTemporal*2 ;
    fHighNorm = 51/fsTemporal*2 ;
    filterOrder = 2 ;
    
    % filter each channels
    for ichannel = setdiff(1:size(reshapeLFPs,1), badChannels)
        [b,a] = butter(filterOrder,[fLowNorm fHighNorm],'stop');
        filteredLFPs(ichannel,:) = filter(b,a,reshapeLFPs(ichannel,:)) ;
    end
    dataOut = reshape(filteredLFPs,10,10,[]) ;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid = vector2grid(vector)
% Converts a 1x100, 100x1 or 100xt vector into a 10x10 or 10x10xt grid,
% rotated correctly. channelDim optionally specifies the dimension
% specifying channel number, which must have length 100.
%
% NOTE: THIS TRANSFORMATION CHANGES THE INDEXING OF THE DATA
% A row vector is stored as [1 2 3 4 ... 99 100]. However, in a 10x10 grid
% format, Paul's group stores the configuration of electrodes as:
%       [100  ......  10 ]
%       | .            . |
%       | .            . |
%       | .            2 |
%       [91   ......   1 ]
% This is 180 degrees rotated from the default MATLAB reshaping. The
% relationship between the index of a value is a row or column vector IV
% and the index of the same value in this rotated matrix IM is:
%    IV = 101 - IM
sv = size(vector);
nrow = sv(1);
ncol = sv(2);
ngrid = 10;

% Check which dimension has length 100 (preferentially choosing the first
% dimension)
if nrow == 100
    t = sv(end);
elseif ncol == 100
    t = nrow;
    vector = vector';
elseif nrow == 10 && ncol == 10
    return
elseif sqrt(nrow) == floor(sqrt(nrow))
    % NEW FUNCTIONALITY: If none of these cases apply, use the first
    % dimension as indicating channel number, as long at is a square
    t = ncol;
    ngrid = sqrt(nrow);
else
    error('vectorLength', 'Invalid input vector size!')
end

% Reshape to 10x10 grid and rotate
if t==1
    grid = rot90(reshape(vector,ngrid,ngrid), 2);
else
    % Reshape to 10x10 grid at every time step
    grid = zeros([ngrid,ngrid,sv(2:end)]);
    for it = 1:prod(sv(2:end))
        grid(:,:,it) = rot90(reshape(vector(:,it),ngrid,ngrid), 2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%