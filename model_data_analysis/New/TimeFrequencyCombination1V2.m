function TimeFrequencyCombination1V2(R,seg_given)
% adjust from function TimeFrequencyCombination1.m and work on mean LFP
% signal
% 1: theta rhythm
% 2: LFP signal
% 3: wavelet spectrogram
[R] = GetBurst2(R);
dt = R.dt;
step_tot = R.step_tot;
seg_size = 2e4; % 1s
seg_num = ceil(step_tot/seg_size);
frequency = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;
pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies
AnalyticP = R.LFP.LFP_gamma_hilbert_abs.^2;
tic;
if nargin == 2
    seg_num = 1;
end
for seg = 1:seg_num
    seg_ind = get_seg(step_tot, seg_size, seg);
    if nargin == 3
        seg_ind = seg_given;
    end
    t = (seg_ind-1)*dt*1e-3; % second
    
    dt = R.dt;
    fs = 1/(dt*1e-3); % sampling frequency (Hz)
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 3; % theta band
    hiFreq = 4;
    Wn = [lowFreq hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
    LFP = mean(R.LFP.LFP{1});
    R.LFP.LFP_theta = filter(b,a,LFP);
    
    h = subplot(3,1,1);
    window_size = 500; % 5ms
    LFP_theta = tsmovavg(R.LFP.LFP_theta,'s',window_size,2);
    plot(t,LFP_theta(seg_ind))
    ylabel('LFP theta rhythm')
    %         [~,locs1] = findpeaks(LFP_theta(i,seg_ind));
    %         [~,locs2] = findpeaks(-LFP_theta(i,seg_ind));
    %         l1 = length(locs1);
    %         l2 = length(locs2);
    p = get(h,'position');
    %         if seg == 10
    %             a = 0;
    %         else
    %             a = 1;
    %         end
    %         raster_plot_YL(R,1,1,[],ceil(dt*seg_ind(1)):ceil(dt*seg_ind(end))+a)
    %         text(-0.15,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    %         ys = zeros(size(x));
    %         ye = 500*ones(size(x));
    %         for j = 1:length(x)
    %             hold on;
    %             plot([x(j) x(j)],[ys(j) ye(j)],'--r')
    %         end
    
    x_tmp = mean(R.LFP.LFP_gamma(:,seg_ind));
    coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
    CData = transpose(coeffs_tmp); % coeffs_tmp'
    subplot(3,1,2) % (2,1,1)
    %         for j = 1:no
    %             plot(t,log10(AnalyticP(j,seg_ind)));
    %             hold on;
    %         end
    %         ylabel('loh_{10}Analytic Power')
    plot(t,mean(R.LFP.LFP_broad(:,seg_ind)),'color',[0.8 0.8 0.8])
    hold on;
    plot(t,x_tmp,'k')
    legend('broadband','gamma-band')
    %         text(-0.1,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    ylabel('LFP(a.u.)') % 14*LFP ('LFP(uV)')
    
    subplot(3,1,3) % (2,1,2)
    uimagesc(t,pseudoFreq(end:-1:1),CData(end:-1:1,:));
    colorbar('off');
    xlim([2*(seg-1) 2*seg]);
    ylim(frequency);
    xlabel('Time(s)');
    ylabel('Frequency(Hz)');
    
    hold on;
    ax=axes('Position',[p(1) 0 p(3) 1],'Unit','normalize',...
        'parent',1);
    %         plot(ax,[dt*1e-3*locs1(1:end);dt*1e-3*locs1(1:end)],[-3*ones(1,l1);15*ones(1,l1)],'g--');
    %         hold on;
    %         plot(ax,[dt*1e-3*locs2(1:end);dt*1e-3*locs2(1:end)],[-3*ones(1,l2);15*ones(1,l2)],'r--');
    start = R.LFP.GammaBurstEvent.burst_start_steps;
    ending = start-1 + R.LFP.GammaBurstEvent.burst_du_steps;
    l = length(start);
    plot(ax,[dt*1e-3*start;dt*1e-3*start],[-3*ones(1,l);15*ones(1,l)],'r--');
    hold on;
    plot(ax,[dt*1e-3*ending;dt*1e-3*ending],[-3*ones(1,l);15*ones(1,l)],'g--');
    xlim([2*(seg-1) 2*seg])
    set(ax,'Xtick',[],'Ytick',[],'Visible','off');
    
    next = input('\t Next figure?');
    delete(gcf);
end
toc;
end