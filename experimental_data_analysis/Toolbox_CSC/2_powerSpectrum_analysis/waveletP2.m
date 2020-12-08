function waveletP2(sigOri,fsTemporal)
% generate wavelet spectrogram for project 2

%% bursting phenomena in wavelet
% for lower frequence bands(<20Hz)
figure;
timeStart = 40 ;
timeEnd = 60 ;
tempData = squeeze (sigOri(5,5,timeStart*fsTemporal+1:timeEnd*fsTemporal)) ;
[wt1,f12,coi] = cwt(tempData,fsTemporal,'VoicesPerOctave',30);
tempTimeAxis = linspace(timeStart,timeEnd, size(tempData,1));
uimagesc(linspace(0,20, size(tempData,1)),flip(f12(135:end)),(zscore(abs(wt1(end:-1:135,:).^2),[],2)))
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar

%% for higher frequence bands(<20Hz)
figure;
timeStart = 8 ;
timeEnd = 10 ;
tempData = squeeze (sigOri(5,5,timeStart*fsTemporal+1:timeEnd*fsTemporal)) ;
[wt2,f22,coi] = cwt(tempData,fsTemporal,'VoicesPerOctave',30);
tempTimeAxis = linspace(timeStart,timeEnd, size(tempData,1));
uimagesc(linspace(0,2, size(tempData,1)),flip(f22(35:end)),(zscore(abs(wt2(end:-1:35,:).^2),[],2)))
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar
end