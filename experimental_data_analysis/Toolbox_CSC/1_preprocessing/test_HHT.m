%% Hilbert Huang Transform

%% preprocess the raw LFP data
fsTemporal = Fs ;
flagBandstop = 0 ;
[sigOri,~,badChannels] = preprocess_LFP(LFPs, flagBandstop) ;

%% 
close all
X = squeeze(sigOri(5,5,1*fsTemporal:10*fsTemporal)) ;
[imf] = emd(X);
%%
Ts = 1/fsTemporal ;
b = [] ;
d = [] ;
for k = 1:length(imf)
   b(k) = sum(imf{k}.*imf{k});
   th   = angle(hilbert(imf{k}));
   d{k} = diff(th)/Ts/(2*pi);
end
[u,v] = sort(-b);
b     = 1-b/max(b);

% Set time-frequency plots.
N = length(X);
c = linspace(0,(N-2)*Ts,N-1);
for k = v(1:end)
   figure, plot(c,d{k},'k.')%,'Color',b([k k k]),'MarkerSize',3);
   set(gca,'FontSize',8,'XLim',[0 c(end)],'YLim',[0 1/2/Ts]); xlabel('Time'), ylabel('Frequency');
end

%% Set IMF plots.
M = length(imf);
N = length(X);
c = linspace(0,(N-1)*Ts,N);
for k1 = 0:4:M-1
   figure
   for k2 = 1:min(4,M-k1), subplot(4,1,k2), plot(c,imf{k1+k2}); set(gca,'FontSize',8,'XLim',[0 c(end)]); end
   xlabel('Time');
end