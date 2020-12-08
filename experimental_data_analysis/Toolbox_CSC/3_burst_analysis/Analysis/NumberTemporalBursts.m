% mean number of temporal bursts

figure;
selectM = [0,1,10] ;
for i = 1:3
    surMethodNum = selectM(i) ;
%% Generate surrogate data (optional)
surSig = generateSur(sigOri,surMethodNum,badChannels) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find subband signals
% apply bandpass filtering for Gamma signals
subBand = [30,80] ;
bandpassSig = find_bandpassSig(surSig,subBand, fsTemporal,3,badChannels) ;
% Hilbert transform for analytic signals
hilbertSig = find_Hilbert(bandpassSig, fsTemporal,4) ;
% find amplitdue of the analytic Gamma as the input
sigIn = abs(squeeze((hilbertSig))) ;

%% Find Gamma burst in the time domain for fig.1
numChannels = size(sigIn,1)*size(sigIn,2) ;
sigReshape = reshape(sigIn,numChannels,[]) ;
stdVal = 3;
GammaBurstEvent = find_Burst_1D(sigReshape,fsTemporal,0,badChannels,...
    numChannels,stdVal) ;

size(GammaBurstEvent.burst_du_steps{25})

sumSize = 0 ;
for i2 = 1:99
    sumSize = sumSize+length(GammaBurstEvent.burst_du_steps{i2}) ;
end
sumSize/99

[x,n] = histcounts(GammaBurstEvent.burst_du_steps{55}) ; loglog(n(1:end-1),x,'o')
hold on
end