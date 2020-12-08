function [dataOut,fsTemporal, badChannels] = preprocess_wholeBrain_1(...
    dataIn, meta, flagZscore)

% [dataOut,fsTemporal, badChannels] = preprocess_wholeBrain_1(...
%    dataIn, meta, flagZscore)
%
% This function implements the procdure to load the anethesia whole brain 
% data from Thomas's group. 
% 
% Inputs:
%       dataIn       (ChannelX X ChannelY X Time)       
%       meta         a structure containing mask information
%       flagZscore   data normalisation, default = 1;
%
% Outputs:
%       dataOut      (ChannelX X ChannelY X Time)
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default variables
if ~exist('flagZscore','var')
    flagZscore = 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale by which to spatially downsample data
mouseResize = 1;
% Flip data upside down, and keep only the right hemisphere
if flagZscore
    data = zscore(rot90(dataIn(20:115, 81:144, :), 2));
else
    data = rot90( dataIn(20:115, 81:144, :),2 ) ;
end
% Spatially downsample each data frame
firstFrame = imresize(data(:,:,1), mouseResize);
% Do the same to the supplied mask of the right cortical hemisphere
cortexMask = imresize(rot90(meta.traceMask(20:115, 81:144,:),2), ...
    mouseResize);
% Output matrix initialisation
dataOut = zeros([size(firstFrame), size(data, 3)]);
% Extract the cortex data
for itime = 1:size(data,3)
    dataOut(:,:,itime) = imresize(data(:,:,itime), mouseResize).*cortexMask;
end
% Sampling rate in time
fsTemporal = 50;
badChannels = find(cortexMask == 0) ;