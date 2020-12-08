function TimeFrequencyCombination1(R,IndC,LFP_centre_x,LFP_centre_y,no_given,seg_given)
% combination of time-LFP & TimeFrequencyCImage.m
% 1: raster plot
% 2: theta/beta/gamma power
% 3: wavelet spectrogram
[R] = GetBurst(R);
dt = R.dt;
step_tot = R.step_tot;
seg_size = 1e4; % 2s
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_gamma);
% frequency = [R.LFP.lowFreq R.LFP.hiFreq];
% scales = R.LFP.wavelet.scales;
% pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies 80 ~ 30 Hz
nos = 1:no;
if nargin >= 5
    nos = no_given;
end
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
% AnalyticP = R.LFP.LFP_gamma_hilbert_abs.^2;
tic;
for i = nos
    if nargin == 6
        seg_num = 1;
    end
    for seg = 1:seg_num
        seg_ind = get_seg(step_tot, seg_size, seg);
        if nargin == 6
            seg_ind = seg_given;
        end
        t = (seg_ind-1)*dt*1e-3; % second
        fs = 1/(dt*1e-3); % sampling frequency (Hz)
        
        h = subplot(2,1,1);
        if seg == seg_num
            a = 0;
        else
            a = 1;
        end
        RasterPlotYL(R,1,1,[],ceil(dt*seg_ind(1)):ceil(dt*seg_ind(end))+a)
        %         text(-0.1,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        subplot(2,1,2)
        % Butterworth filter
        order = 4; % 4th order
        lowFreq = [4 20 30]; % theta/beta/gamma
        hiFreq = [7 30 80];
        LFP = R.LFP.LFP{1}(:,seg_ind);
        ca = cell(1,length(IndC)*3);
        for k = 1:3
            Wn = [lowFreq(k) hiFreq(k)]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            for j = 1:length(IndC)
                [~,indE] = sort(Distance_xy(LFP_centre_x,LFP_centre_y,Lattice(IndC(j),1),Lattice(IndC(j),2),2*hw));
                LFP_tbg = filter(b,a,LFP(indE(1),:));
                LFP_power = abs(hilbert(LFP_tbg)).^2/(hiFreq(k)-lowFreq(k));
                ca{(k-1)*length(IndC)+j} = sprintf('Electrode No.%d',indE(1));
%                 plot(t,LFP_power);
                switch k
                    case 1
                        plot(t,LFP_power,'g','LineWidth',0.75*j); % 'Color',[0+0.2*(j-1) 1 0]) % theta green
                    case 2
                        plot(t,LFP_power,'b','LineWidth',0.75*j); % 'Color',[0 0 1-0.2*(j-1)]) % beta blue
                    case 3
                        plot(t,LFP_power,'r','LineWidth',0.75*j); % 'Color',[1 0 0+0.2*(j-1)]) % gamma red 
                end
                hold on;
            end
        end
        legend(ca);
        xlabel('Time(s)')
        ylabel('Power(mV^2/Hz)')
        %         window_size = 500; % 5ms
        %         LFP_theta = tsmovavg(LFP_theta,'s',window_size,2);
        %         plot(t,LFP_theta(i,seg_ind))
        %         ylabel('LFP theta rhythm')
        %         [~,locs1] = findpeaks(LFP_theta(i,seg_ind));
        %         [~,locs2] = findpeaks(-LFP_theta(i,seg_ind));
        %         l1 = length(locs1);
        %         l2 = length(locs2);
%         ind = input('\t Please select electrode index:');
        p = get(h,'position');
%         % Butterworth filter
%         order = 4; % 4th order
%         lowFreq_br = 5; % broad band (1-1000 Hz)
%         hiFreq_br = 100;
%         Wn = [lowFreq_br hiFreq_br]/(fs/2);
%         [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
%         LFP_tbg = filter(b,a,R.LFP.LFP{1}(ind,:));
%         x_tmp = LFP_tbg(seg_ind);
%         %         coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
%         fc = centfrq('cmor1.5-1');
%         scalerange = fc./([5 100]*(1/fs));
%         scales = scalerange(end):0.5:scalerange(1);
%         pseudoFreq = scal2frq(scales,'cmor1.5-1',1/fs);
%         coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % this had narrow width of scales
%         CData = transpose(coeffs_tmp); % coeffs_tmp'
%         subplot(3,1,3)
%         %         for j = 1:no
%         %             plot(t,log10(AnalyticP(j,seg_ind)));
%         %             hold on;
%         %         end
%         %         ylabel('loh_{10}Analytic Power')
%         %         plot(t,R.LFP.LFP_broad(i,seg_ind),'color',[0.8 0.8 0.8])
%         %         hold on;
%         %         plot(t,x_tmp,'k')
%         %         legend('broadband','gamma-band')
%         %         text(-0.1,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%         %         ylabel('LFP(a.u.)') % 14*LFP ('LFP(uV)')
%         %
%         %         subplot(4,1,4) % (3,1,3)
%         %         zscoreCoef = zscore(CData,1);
%         %         uimagesc(t,pseudoFreq(end:-1:1),zscoreCoef(end:-1:1,:));
%         uimagesc(t,pseudoFreq(end:-1:1),CData(end:-1:1,:));
%         colorbar('off');
%         xlim([1*(seg-1) 1*seg]);
%         ylim([5 100]);
%         xlabel('Time(s)');
%         ylabel('Frequency(Hz)');
%         set(gca,'YDir','normal')
        
        %         subplot(4,1,4)
        %                 zscoreCoef = CData;
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
        %         temp = zeros(size(binaryImage(:,1:2*fsTemporal)));
        %         for iRegion = 1:length(validRegionIdx)
        %             temp(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
        %         end
        %         uimagesc(t,pseudoFreq(end:-1:1),temp(end:-1:1,:))
        %         % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
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
        %         start = R.LFP.GammaBurstEvent.burst_start_steps{i};
        %         ending = start-1 + R.LFP.GammaBurstEvent.burst_du_steps{i};
        %         l = length(start);
        %         plot(ax,[dt*1e-3*start;dt*1e-3*start],[-3*ones(1,l);15*ones(1,l)],'r--');
        %         hold on;
        %         plot(ax,[dt*1e-3*ending;dt*1e-3*ending],[-3*ones(1,l);15*ones(1,l)],'g--');
        %         xlim([0 2]) % ([2*(seg-1) 2*seg])
        set(ax,'Xtick',[],'Ytick',[],'Visible','off');
        next = input('\t Next figure?');
        delete(gcf);
    end
    toc;
end
end
