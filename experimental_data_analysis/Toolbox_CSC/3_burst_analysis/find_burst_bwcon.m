function [CC,validRegionIdx, imageBurst, centroids, boundary] = ...
    find_burst_bwcon(sigIn, prctileDetect,minConnect)

% [CC,validRegionIdx, imageBurst, centroids, boundary] = ...
%    find_burst_bwcon(sigIn, prctileDetect,minConnect)
%
% This function implements ...
% 
% Input:
%       sigIn            multi-dimensional signals
%       prctileDetect    percentiles for threshold
%       minConnect       minimum connection required for valid regions
%
% Output:
%       CC               (ChannelX X ChannelY X Time) matrix of 
%       validRegionIdx   Index for valid regions
%       imageBurst       new binary image with valid regions
%       centroids        centroids of valid regions
%       boundary         boundaries of valid regions
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings
if ~exist('pencentile','var')
    prctileDetect = 95 ;  
end
if ~exist('minConnect','var')
    minConnect = 200 ;   
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find burst 2D
Y = prctile(zscoreCoef(:),prctileDetect) ;
binaryImage = zeros(size(sigIn)) ;
binaryImage(sigIn>=Y) = 1 ;

%% find connected regions
CC = bwconncomp(binaryImage(:,:)) ;
validRegionIdx = 0;
count = 1;
for iRegion = 1:size(CC.PixelIdxList,2)
    if(size(CC.PixelIdxList{iRegion},1)>minConnect)
        validRegionIdx(count) = iRegion ;
        count = count + 1 ;
    end
end

%% new images with valid regions
imageBurst = binaryImage ;
for iRegion = 1:length(validRegionIdx)
    imageBurst(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
end


%% central points of valid regions
S = regionprops(CC,'Centroid');
centroids = cat(1, S.Centroid);

%% boundary of valid regions
B = regionprops(CC,'BoundingBox');
boundary = cat(1, B.BoundingBox);

