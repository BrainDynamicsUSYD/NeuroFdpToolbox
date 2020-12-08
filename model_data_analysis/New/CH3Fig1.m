figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'CH3Fig1', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

disp('plot_SWR...');

dt = R.dt;
step_tot = R.step_tot;
seg_size = 1*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_broad);
nos = 1:no;
for i = nos
    for seg = 4:seg_num
        seg_ind = get_seg(step_tot, seg_size, seg);
        t = (seg_ind-1)*dt*1e-3; % second
        ax1 = subplot(4,1,1);
        LFP_tmp = R.LFP.LFP_broad(i,seg_ind);
        plot(t, LFP_tmp);
        ylabel('Broad LFP','fontsize',10)
        text(-0.1,1,'A','Units', 'Normalized','FontSize',12)
        ax2 = subplot(4,1,2);
        rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
        plot(t,  rip_tmp, 'b');
        ylabel('Ripple LFP','fontsize',10);
        hold on;
        hil_tmp = R.LFP.LFP_ripple_hilbert(i,seg_ind);
        plot(t, hil_tmp,'r'); 
%         plot(t, ones(size(t))*R.LFP.ripple_event.hil_mean_baseline(i,end),'m--');
%         plot(t, ones(size(t))*(R.LFP.ripple_event.hil_mean_baseline(i,end) + R.LFP.ripple_event.no_std*std(R.LFP.LFP_ripple_hilbert(i,:))),'m-');
%         plot(t, ones(size(t))*(R.LFP.ripple_event.hil_mean_baseline(i,end) + R.LFP.ripple_event.peak_no_std*std(R.LFP.LFP_ripple_hilbert(i,:))),'m-.');
        scales = R.LFP.wavelet.scales;
        x_tmp = R.LFP.LFP_ripple(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
        ax3 = subplot(4,1,3); % Scaleogram with pseudo-Frquency
        freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
        imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
        ylim(freqrange)
        ylabel('Frequency(Hz)','fontsize',10)
        text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
        ax4 = subplot(4,1,4);
        reduced = R.reduced;
        s_tmp = R.ExplVar.LFP_range_sigma;
        spike_sort_range = 10;
        spike_sort_neurons = R.LFP.LFP_neurons{1}(i,:) >= 1/(s_tmp*sqrt(2*pi))*exp(-0.5*(spike_sort_range/s_tmp)^2);
        reduced.spike_hist{1} = reduced.spike_hist{1}(spike_sort_neurons, :);
        R_LFP.N(1) = sum(spike_sort_neurons);
        R_LFP.reduced = reduced; clear reduced;
        dt_conv = R.dt/R.reduced.dt;
        raster_plot(R_LFP, 1, seg, [], 'seg_size', round(seg_size* dt_conv))
        xlabel('time(s)','fontsize',10);
        text(-0.1,1,'D','Units', 'Normalized','FontSize',12)
        % highlight detected ripple events
        if ~isempty(R.LFP.ripple_event.ripple_du_steps{i})
            ripple_start = R.LFP.ripple_event.ripple_start_steps{i};
            ripple_du = R.LFP.ripple_event.ripple_du_steps{i};
            hold on;
            y_lim = get(gca,'ylim');
            for i_r = 1:length(ripple_du)
                tA = ripple_start(i_r);
                tB = tA + ripple_du(i_r) - 1;
                plot([tA*dt tA*dt]*1e-3,y_lim,'r');
                plot([tB*dt tB*dt]*1e-3,y_lim,'r');
            end
        end
        % Link axes to synchronise them when zooming
        linkaxes([ax1 ax2 ax3 ax4],'x');
        set(ax1, 'xlim', minmax(t));
        next = input('\t Next figure?');
        delete(gcf);
    end
end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig1 % this is the trick!!