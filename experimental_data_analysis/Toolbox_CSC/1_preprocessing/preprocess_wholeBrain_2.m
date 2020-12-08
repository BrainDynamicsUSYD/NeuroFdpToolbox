function [dataOut,fsTemporal, badChannels] = preprocess_wholeBrain_2...
    (dataIn, mask, flagZscore, flagHighpass, flagGaussianFilt,assignZero)

% [dataOut,fsTemporal, badChannels] = preprocess_wholeBrain_2...
%    (dataIn, mask, flagZscore, flagHighpass, flagGaussianFilt)
%
% This function implements the procdure to load the awkae type whole brain 
% data from Thomas's group and includes a highpass filter and a spatial 
% gaussian filter. 
% 
% Inputs:
%       dataIn      (ChannelX X ChannelY X Time)
%       mask        (ChannelX X ChannelY X Time) of Boolean for brain areas
%       flags       default = 1;
%
% Outputs:
%       dataOut     (ChannelX X ChannelY X Time)
%       fsTemporal
%       badChannels
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default variables
if ~exist('flagZscore','var')
    flagZscore = 1;
end
if ~exist('flagHighpass','var')
    flagHighpass = 1;
end
if ~exist('flagGaussianFilt','var')
    flagGaussianFilt = 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zscore the data
reshapeIn = reshape(dataIn,size(dataIn,1)*size(dataIn,2),[]) ;
for iChannel = 1:size(reshapeIn,1)
    if flagZscore==1
        reshapeIn = zscore(reshapeIn(iChannel,:)) ;
    elseif flagZscore==2
        reshapeIn = reshapeIn(iChannel,:) - mean(reshapeIn(iChannel,:)) ;
    end
end
dataOut = zeros(size(dataIn)) ;

% apply high pass filter (frequency larger than 0.1Hz)
if flagHighpass
    [B,A] = cheby2(5,20,0.1/150*2,'high') ;
    [sizeMaskX, sizeMaskY] = size(mask) ;
    
    % loop through all the channels
    for xNode = 1:sizeMaskX
        for yNode = 1:sizeMaskY
            if mask(xNode,yNode) == 1
                dataOut(xNode,yNode,:) = filtfilt( B,A,double( squeeze( ...
                    dataIn(xNode,yNode,:)) ) );
            else
                if assignZero == 1
                    dataOut(xNode,yNode,:) = zeros(size(dataIn,3),1);
                else
                    dataOut(xNode,yNode,:) = nan(size(dataIn,3),1);
                end
            end
        end
    end
else
    [sizeMaskX, sizeMaskY] = size(mask) ;
    % loop through all the channels
    for xNode = 1:sizeMaskX
        for yNode = 1:sizeMaskY
            if mask(xNode,yNode) == 1
                dataOut(xNode,yNode,:) = dataIn(xNode,yNode,:) ;
            else
                if assignZero == 1
                    dataOut(xNode,yNode,:) = zeros(size(dataIn,3),1);
                else
                    dataOut(xNode,yNode,:) = nan(size(dataIn,3),1);
                end
            end
        end
    end
end

if flagGaussianFilt
    %apply spatial smoothing (Gaussian filter)
    for itime = 1:size(dataOut,3)
        dataOut(:,:,itime) = imgaussfilt(dataOut(:,:,itime),0.8).*mask;
    end
end

% Sampling rate in time
fsTemporal = 150 ;
badChannels = find(mask == 0);