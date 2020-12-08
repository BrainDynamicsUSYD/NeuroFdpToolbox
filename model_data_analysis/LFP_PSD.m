function LFP_PSD(R)
% plot average LFP PSD(fft) on log scale
Fs = 1e4; %
% LFP = mean(R.LFP.LFP{1}(:,Fs+1:end)); % discard 1s transient time
LFP = mean(R.LFP.LFP_broad(:,Fs+1:end));
% LFP = R.LFP.LFP{1}(11,Fs+1:end); % discard 1s transient time
% LFP = vec2mat(R.num_spikes{1},10); % discard 1s transient time
% LFP = sum(LFP,2);
% LFP = LFP'/3969*1e3; % Hz on single neuron

%%% FFT estimation %%%
N = length(LFP);
xdft = fft(LFP);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;

% a = 2.5; % 2.96
% y = 1e3./(freq.^a);
window_size = 10; % 1 Hz
psdx = tsmovavg(psdx,'s',window_size,2);
% subplot(1,2,1)
% plot(freq(51:end),psdx(51:end))
% title('PSD of Population Firing Rate')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% subplot(1,2,2)
loglog(freq,psdx)
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
%%
% plot(freq(12:round(end/20)),psdx(12:round(end/20))) % 12:2e3
ind = find(freq > 300);
% freqlog = log10(freq);
% psdxlog = log10(psdx);
f = freq(ind);
z = psdx(ind);
OneOverF = fittype(@(coe,slope,f) coe./(f.^slope),'independent','f','dependent','z');
[Sfit,~] = fit(f(:),z(:),OneOverF,'StartPoint',[1e-4 1.5],'Lower',[1e-4 0.1],'Upper',[10 3]);
Sfit.coe
% Sfit.slope % = 2;
y = (Sfit.coe)./(freq.^(Sfit.slope));
hold on
loglog(freq,y)
% xlim([1e-1 1e3])
% text(10,1e-3,'slope -2.96')
str = ['slope = ',num2str(Sfit.slope)];
text(40,max(psdx),str)
title('PSD of LFP using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
% text(-0.1,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
end