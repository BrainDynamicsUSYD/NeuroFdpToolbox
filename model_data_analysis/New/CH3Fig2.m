figure_width = 11.4; % cm
figure_hight = 9; % cm
figure('NumberTitle','off','name', 'CH3Fig2', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

subplot(2,1,1)
for i = 1:2
    switch i
        case 1
            Fs = 1e4;
            Signal = mean(R.LFP.LFP_broad);
        case 2
            Fs = 1e3;
            Signal = R.reduced.num_spikes{1}/R.N(2)*Fs;
    end    
    %%% FFT estimation %%%
    N = length(Signal);
    xdft = fft(Signal);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/N:Fs/2;
    if i == 1
        window_size = 200; % 5 Hz
        psdx = tsmovavg(psdx,'s',window_size,2);
        yyaxis left
        plot(freq,psdx)        
        ylabel('LFP PSD','fontsize',10)
    else
        window_size = 200; % 5 Hz
        psdx = tsmovavg(psdx,'s',window_size,2);
        yyaxis right
        plot(freq,psdx)
        ylabel('Firing rate PSD','fontsize',10)
    end
    xlabel('Frequency(Hz)','fontsize',10)
    xlim([50 300])
end
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)

sw_amp = cell2mat( R.LFP.wavelet.peak.sw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.sw_amp)) );
rp_freq = cell2mat( R.LFP.wavelet.peak.rp_freq(~cellfun(@isempty,R.LFP.wavelet.peak.rp_freq)) );
rp_raw_amp = cell2mat( R.LFP.wavelet.peak.rp_raw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.rp_raw_amp)) );
dat1 = [sw_amp(~isnan(sw_amp))', rp_raw_amp(~isnan(sw_amp))'];
dat2 = [sw_amp(~isnan(sw_amp))', rp_freq(~isnan(sw_amp))'];
subplot(2,2,3)
n = hist3(dat1,[30 30]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat1(:,1)),max(dat1(:,1)),size(n,1)+1);
yb = linspace(min(dat1(:,2)),max(dat1(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
[r,p] = corrcoef(sw_amp(~isnan(sw_amp)), rp_raw_amp(~isnan(sw_amp)));
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([6 176 176 6 6],[11 11 35 35 11],'k')
xlabel('Wave Amplitude(a.u.)','fontsize',10)
ylabel('Ripple Amplitude(a.u.)','fontsize',10)
xlim([6 176])
ylim([11 35])
text(-0.28,1,'B','Units', 'Normalized','FontSize',12)

subplot(2,2,4)
n = hist3(dat2,[30 30]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat2(:,1)),max(dat2(:,1)),size(n,1)+1);
yb = linspace(min(dat2(:,2)),max(dat2(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
[r,p] = corrcoef(sw_amp(~isnan(sw_amp)), rp_freq(~isnan(sw_amp)));
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([6 176 176 6 6],[100 100 166 166 100],'k')
xlabel('Wave amplitude(a.u)','fontsize',10)
ylabel('Ripple Frequency(Hz)','fontsize',10)
xlim([6 176])
ylim([100 166])
text(-0.28,1,'C','Units', 'Normalized','FontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig2 % this is the trick!!