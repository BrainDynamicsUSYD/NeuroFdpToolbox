%% This function finds the peak amplitude and intervals of Gamma oscillation
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
tempData = squeeze(sigOri(5,5,:)) ;
tempData = rand(10000,1) ;
% tempData = squeeze(dataFull(5,5,:)) ;
N = 4;  % Filter order
% discardTimeSecs = 5;
% Filter LFPs
disp('Filtering waveforms...') ; tic
fLow = 20 ;
fHigh = 80 ;
% Design Butterworth band-pass filter
h = fdesign.bandpass('N,F3dB1,F3dB2',N,fLow,fHigh,fsTemporal);
Hd = design(h, 'butter');
set(Hd, 'Arithmetic', 'double');

SOS = Hd.sosMatrix;
G = Hd.ScaleValues;

 % Filter signals forwards and backwards to avoid phase distortion
bandpassLFPs = filtfilt(SOS,G,tempData)';        



ICI = [] ;
Amp = [] ;
tempAmp = [];
    tempICI = [];
    [p1,l1] = findpeaks(bandpassLFPs);
    [p2,l2] = findpeaks(-bandpassLFPs);
    l = min(length(l1),length(l2));
    m = 1; % index for l2
    for j = 1:l
        while m <= length(l2) && l1(j) >= l2(m)
            m = m + 1;
        end
        if m > length(l2)
            break
        end
        tempAmp = [tempAmp p1(j) + p2(m)];
        m = m + 1;
    end
    l = length(tempAmp);
    tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
    if length(tempAmp) < length(tempICI)
        disp('Miss Match!');
    end
%
dir = 1 % forward; 2 for backward
switch dir
    case 1
        tempAmp = tempAmp(1:length(tempICI));
    case 2
     tempAmp = tempAmp(2:end);
     tempICI = tempICI(1:length(tempAmp));
end
    ICI = [ICI tempICI];
    Amp = [Amp tempAmp];
    
dat = [Amp', ICI'];
%%
figure
nFFT = hist3(dat,[32 32]);
    n1 = nFFT';
    n1(size(nFFT,1) + 1, size(nFFT,2) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(nFFT,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(nFFT,1)+1);
    h = pcolor(xb,yb,n1);
    set(h, 'EdgeColor', 'none');
    h.ZData = ones(size(n1)) * -max(max(nFFT));
    colormap(hot)
    oldcmap = colormap(gray);
    colormap( flipud(oldcmap) );
    xlabel('amplitude')
    ylabel('ICI (ms)')
   
    [r,p] = corrcoef(Amp',ICI);
 
if p(2) > 10^(-4)
    pvalue = ['p = ',num2str(p(2),'%.2f')];
else
    pvalue = ['p < 10^{-4}'];
end
%     title(['r = ',num2str(r(2),'%.2f'),pvalue])
tp = {['r = ',num2str(r(2),'%.2f')],pvalue}
