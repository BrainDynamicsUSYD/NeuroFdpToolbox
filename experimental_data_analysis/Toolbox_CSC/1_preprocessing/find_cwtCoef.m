function [wvcfs,cfreqso] = find_cwtCoef(sigIn,cfreq, fsTemporal,morletParam)

% [wvcfs,cfreqso] = find_cwtCoef(sigIn,cfreq, fsTemporal) ;
%
% This function implements finding the complex wavelet coefficient by
% calling Rory's code (cwtft).
% 
% Input:
% sigIn          1x1xT vector time series or multi-dimensional vector with 
%                   time specified in dimension DIM
% fsTemporal     sampling frequnecy
% cfreq          centre frequency
%
% Output:
% wvcfs          (Freq X ChannelX X ChannelY X Time) matrix of wavelet
%                coefficient
% cfreqso        
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings

addpath(genpath([pwd,'/ToolNeuroPatt/NeuroPattToolbox-master']))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering parameters
tic
% morletParam = 6 ;
timeDim = 3 ;

% Filter LFP signals
disp('Filtering waveforms...') ; tic
[cfx,cfreqso] = morletWaveletTransform(sigIn, fsTemporal, cfreq, ...
    morletParam, timeDim) ;
wvcfs = squeeze(cfx) ;
toc

rmpath(genpath([pwd,'/ToolNeuroPatt/NeuroPattToolbox-master']))
