function TimeFrequencyCombination2(R,no_given,seg_given)
% combination of time-LFP & TimeFrequencyCImage.m
% 1: theta rhythm
% 2: LFP signal
% 3: wavelet spectrogram
% 4: binary wavelet spectrogram
[R] = GetBurst(R);
dt = R.dt;
step_tot = R.step_tot;
seg_size = 2e4; % 1s
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_gamma);
frequency = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;
pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies 80 ~ 30 Hz
nos = 1:no;
if nargin > 2
    nos = no_given;
end
% AnalyticP = R.LFP.LFP_gamma_hilbert_abs.^2;
tic;
for i = nos
    if nargin == 4
        seg_num = 1;
    end
    for seg = 2:seg_num
        seg_ind = get_seg(step_tot, seg_size, seg);
        if nargin == 4
            seg_ind = seg_given;
        end
        t = (seg_ind-1)*dt*1e-3; % second
        fs = 1/(dt*1e-3); % sampling frequency (Hz)
        
        % Butterworth filter
        %         order = 4; % 4th order
        %         lowFreq = 3; % theta band
        %         hiFreq = 4;
        %         Wn = [lowFreq hiFreq]/(fs/2);
        %         [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
        %         LFP = R.LFP.LFP{1};
        %         for j = 1:no
        %             LFP_theta(j,:) = filter(b,a,LFP(j,:)); %#ok<AGROW>
        %         end
        %         R.LFP.LFP_theta = LFP_theta;
        
        h = subplot(2,1,1);
        %         window_size = 50; % 5ms
        %         LFP_theta = tsmovavg(LFP_theta,'s',window_size,2);
        %         plot(t,LFP_theta(i,seg_ind))
        %         text(-0.1,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        %         ylabel('LFP theta rhythm')
        %         [~,locs1] = findpeaks(LFP_theta(i,seg_ind));
        %         [~,locs2] = findpeaks(-LFP_theta(i,seg_ind));
        %         l1 = length(locs1);
        %         l2 = length(locs2);
        p = get(h,'position');
        x_tmp = R.LFP.LFP_gamma(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
        CData = transpose(coeffs_tmp); % coeffs_tmp'
        %         subplot(4,1,2)
        plot(t,R.LFP.LFP_broad(i,seg_ind),'color',[0.8 0.8 0.8])
        hold on;
        plot(t,x_tmp,'k')
        legend('broadband','gamma-band')
        %         text(-0.1,1,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        ylabel('LFP(a.u.)') % 14*LFP ('LFP(uV)')
        subplot(2,1,2)
        %         [wt,f,coi] = cwt(R.LFP.LFP{1}(14,seg_ind),fs,'VoicesPerOctave',30);
        %         for j = 1:length(coi)
        %             ind = find(f<=coi(j));
        %             wt(ind,j) = NaN;
        %         end
        %         ind = find(f<30 | f >80);
        %         wt(ind,:) = [];
        %         f(ind) = [];
        %         surface(t,f,abs(wt)) % cfs
        %         axis tight
        %         shading flat
        uimagesc(t,pseudoFreq(end:-1:1),CData(end:-1:1,:));
        colorbar('off');
        %         text(-0.1,0.98,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        xlim([2*(seg-1) 2*seg]);
        ylim(frequency);
        xlabel('Time(s)');
        ylabel('Frequency(Hz)');
        set(gca,'YDir','normal')
        %         subplot(4,1,3)
        %         plot(t,S.I_AMPA(4,seg_ind)+S.I_ext(4,seg_ind))
        %         xlabel('Time(s)')
        %         ylabel('I_{excitatory}(nA)')
        %         subplot(4,1,4)
        %         plot(t,S.I_GABA(4,seg_ind))
        %         xlabel('Time(s)')
        %         ylabel('I_{inhibitory}(nA)')
        %         zscoreCoef = CData;
        %         Y = prctile(zscoreCoef(:),95) ;
        %         binaryImage = zeros(size(zscoreCoef)) ;
        %         binaryImage(zscoreCoef>=Y) = 1 ;
        %         CC = bwconncomp(binaryImage(:,1:end)) ;
        %         validRegionIdx = 0;
        %         count = 1;
        %         for iRegion = 1:size(CC.PixelIdxList,2)
        %             if(size(CC.PixelIdxList{iRegion},1)>200)
        %                 validRegionIdx(count) = iRegion ;
        %                 count = count + 1 ;
        %             end
        %         end
        %         fsTemporal = 1e4;
        %         temp = zeros(size(binaryImage(:,1:1*fsTemporal)));
        %         for iRegion = 1:length(validRegionIdx)
        %             temp(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
        %         end
        %         uimagesc(t,pseudoFreq(end:-1:1),temp(end:-1:1,:))
        %         text(-0.1,0.96,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
        %         set(gca,'YDir','normal')
        %         ylabel('Frequency (Hz)')
        %         xlabel('Time (s)')
        %         S = regionprops(CC,'Centroid');
        %         centroids = cat(1, S.Centroid);
        %
        %         cTime = centroids(validRegionIdx(:),1)/fsTemporal ;
        %         cFreq = pseudoFreq(round(centroids(validRegionIdx(:),2))) ;
        %
        %         B = regionprops(CC,'BoundingBox');
        %         boundary = cat(1, B.BoundingBox);
        %
        %         duTime = boundary(validRegionIdx(:),3)/fsTemporal ;
        %         freqLower = min([length(pseudoFreq)*ones(length(validRegionIdx),1),...
        %             round( centroids(validRegionIdx(:),2)+boundary(validRegionIdx(:),4)/2)]') ;
        %         freqUpper = max([1*ones(length(validRegionIdx),1),...
        %             round( centroids(validRegionIdx(:),2)-boundary(validRegionIdx(:),4)/2)]') ;
        %         bwFreq = pseudoFreq(freqUpper ) - pseudoFreq(freqLower) ;
        
        
        hold on;
        ax=axes('Position',[p(1) 0 p(3) 1],'Unit','normalize',...
            'parent',1);
        
        %         plot(ax,[dt*1e-3*locs1(1:end);dt*1e-3*locs1(1:end)],[-3*ones(1,l1);15*ones(1,l1)],'g--');
        %         hold on;
        %         plot(ax,[dt*1e-3*locs2(1:end);dt*1e-3*locs2(1:end)],[-3*ones(1,l2);15*ones(1,l2)],'r--');
        start = R.LFP.GammaBurstEvent.burst_start_steps{i};
        ending = start-1 + R.LFP.GammaBurstEvent.burst_du_steps{i};
        l = length(start);
        plot(ax,[dt*1e-3*start;dt*1e-3*start],[-3*ones(1,l);15*ones(1,l)],'r--');
        hold on;
        plot(ax,[dt*1e-3*ending;dt*1e-3*ending],[-3*ones(1,l);15*ones(1,l)],'g--');
        xlim([2*(seg-1) 2*seg]) % ([0 2]) %
        set(ax,'Xtick',[],'Ytick',[],'Visible','off');
        next = input('\t Next figure?');
        delete(gcf);
    end
    toc;
end
end
