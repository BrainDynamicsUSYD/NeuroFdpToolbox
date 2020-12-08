% load the data
% load the sync info.

%% test evoked activity
binTimeStamps = thisTimeStamp' - fileStartTime ;
binTimeStamps = fix (binTimeStamps*Fs) + 1 ;

%% plot time series of LFP
figure;
testSig = LFPs(20,:) ;
plot(testSig)
hold on
plot(binTimeStamps,max(testSig)*ones(length(binTimeStamps),1),'ro') 

%% fft of the LFP
figure;
freqSig = fft(testSig) ;
lengthSig = length(testSig) ;
freqAxis = linspace(1,Fs,lengthSig) ;

plot(freqAxis(1:lengthSig/2), abs(freqSig(1:lengthSig/2)))
figure;
loglog(freqAxis(1:lengthSig/2), abs(freqSig(1:lengthSig/2)))

%% bandstop electric noise
% bandstop the 49-51 Hz, use 4th order butterworth filter
fLowNorm = 49.75/Fs*2 ;
fHighNorm = 50.25/Fs*2 ;
filterOrder = 2 ;
[b,a] = butter(filterOrder,[fLowNorm fHighNorm],'stop');
filteredLFPs = filtfilt(b,a,testSig) ;

% 
% % bandstop the 49-51 Hz, use 4th order butterworth filter
% fLowNorm = 49.5/Fs*2 ;
% fHighNorm = 50.5/Fs*2 ;
% filterOrder = 2 ;
% [b,a] = butter(filterOrder,[fLowNorm fHighNorm],'stop');
% filteredLFPs = filtfilt(b,a,filteredLFPs) ;

% bandstop the 49-51 Hz, use 4th order butterworth filter
fLowNorm = 99.5/Fs*2 ;
fHighNorm = 100.5/Fs*2 ;
filterOrder = 3 ;
[b,a] = butter(filterOrder,[fLowNorm fHighNorm],'stop');
filteredLFPs = filtfilt(b,a,filteredLFPs) ;

figure;
freqSig = fft(filteredLFPs) ;
lengthSig = length(testSig) ;
freqAxis = linspace(1,Fs,lengthSig) ;

plot(freqAxis(1:lengthSig/2), abs(freqSig(1:lengthSig/2)))

%% try wavelet spectrum
addpath('ToolOthers/uimage')
freqRange = 0.5:100 ;
    fc = centfrq('cmor1.5-1') ;
    scalerange = fc./(freqRange/Fs) ;
    scales = scalerange(end):0.5:scalerange(1) ;
    pseudoFreq = scal2frq(scales, 'cmor1.5-1', 1/Fs) ;
    
    tempData =  squeeze (filteredLFPs(4*Fs:8*Fs)) ; 
    
    wt = cwt( tempData ,scales, 'cmor1.5-1'  ) ;
    tempTimeAxis = linspace(4,8,size(wt,2)) ;
    figure
    wtNorm = zscore(abs(wt),[],2) ;
    
    uimagesc(tempTimeAxis,pseudoFreq(end:-1:1),((abs(wt(end:-1:1,:)) )))
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    set(gca,'YDir','normal')


%% average over stimulus
evokePeriod = [] ;
for evokeIdx = 1: length(binTimeStamps)
    evokePeriod = [evokePeriod,  binTimeStamps(evokeIdx):...
        binTimeStamps(evokeIdx)+fix(0.25*Fs)] ;
end
%%
evokeSig = filteredLFPs(evokePeriod(:)) ;
freqSig = fft(evokeSig) ;
lengthSig = length(evokeSig) ;
freqAxis = linspace(1,Fs,lengthSig) ;

plot(freqAxis(1:lengthSig/2), abs(freqSig(1:lengthSig/2)))
figure;
loglog(freqAxis(1:lengthSig/2), abs(freqSig(1:lengthSig/2)))

%%
% addpath(genpath('ToolOthers/chronux_2_11'))
params.tapers = [4,7] ;
params.Fs = Fs ;
[S,f] = mtspectrumc(evokeSig,params) ;
plot(f,S)
figure;
loglog(f,S)

%% find the evoked activity
[~,startStim] = min(abs(stimTime-0.25)) ;
[~,endStim] = min(abs(stimTime-2)) ;
evokeLFP = mean(arrangedLFP(1,:,startStim:endStim)) ;
params.tapers = [4,7] ;
params.Fs = Fs ;
[S,f] = mtspectrumc(evokeLFP,params) ;
figure;
semilogy(f,S)
figure;
loglog(f,S)

%% find the evoked activity
[~,startStim] = min(abs(stimTime-0.25)) ;
[~,endStim] = min(abs(stimTime-2)) ;
evokeSingle = squeeze(arrangedLFP(1,:,startStim:endStim)) ;
params.tapers = [4,7] ;
params.Fs = Fs ;
S = zeros(100,1025) ;
for iTrial = 1: 100
    [S(iTrial,:),f] = mtspectrumc(evokeSingle(iTrial,:),params) ;
end
S = mean(S) ;
figure;
semilogy(f,S)
figure;
loglog(f,S)


%% find the induced activity
[~,startStim] = min(abs(stimTime-0.25)) ;
[~,endStim] = min(abs(stimTime-2)) ;

inducedLFP = squeeze(arrangedLFP(1,:,startStim:endStim))' - ...
    squeeze(evokeLFP)*ones(1,100) ;
params.tapers = [4,7] ;
params.Fs = Fs ;
% [S,f] = mtspectrumc(inducedLFP(:),params) ;
% figure;
% semilogy(f,S)
% figure;
% loglog(f,S)
%% calculate the average power spectrum
S = zeros(100,1025) ;
for iTrial = 1: 100
    [S(iTrial,:),f] = mtspectrumc(inducedLFP(:,iTrial),params) ;
end
S = mean(S) ;
figure;
semilogy(f,S)
ylim([1e-14 1e-8])
figure;
loglog(f,S)

%% try mtspecgramc
evokeSingle = arrangedLFP(1,1,:) ;
movingwin = [0.4,0.05] ;
params.tapers = [4,7] ;
params.Fs = Fs ;
[s,t,f] = mtspecgramc(evokeSingle,movingwin,params) ;
imagesc(f,t,s)
[~,startStim] = min(abs(t-0.25)) ;
[~,endStim] = min(abs(t-2)) ;
meanSpectrum = mean(s(startStim:endStim,:)) ;
figure;
plot(f,meanSpectrum)





%% try Delta plane wave and Gamma burst




