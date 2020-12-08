figure_width = 11.4; % cm
figure_hight = 5.6; % cm
figure('NumberTitle','off','name', 'CH3Fig7', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

% load('0007-201901141714-05837_in_1547446578141_0_neurosamp.mat', 'ripple_power_grid')
% R = load('0007-201901141714-05837_in_1547446578141_out_RYG.mat');
% R = get_SWR(R);
% S = full(R.spike_hist{1});
% S = reshape(S,63,63,[]);
% 
% mm = minmax(ripple_power_grid(:)');

subplot(1,2,1)
box on;
x = [1 63 63 1];
y = 63-[11 11 21 21];
fill(x,y,'r')
% plot([1,1;63,63],63-[11,21;11,21],'k')
xlim([1 63]);
ylim([1 63]);
set(gca,'xtick',[],'ytick',[]);
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)

subplot(1,2,2)
for t = 31 % 1:length(ripple_power_grid)
    Spikes = sum(S(:,:,(t-25:t+24)+2e4),3);
    [row,col] = find(Spikes);
    imagesc(ripple_power_grid(:,:,t)',mm);
    hold on
    plot([1,1;63,63],[11,21;11,21],'w')
    hold on
    plot(row,col,'k.','MarkerSize',3);    
end
set(gca,'xtick',[],'ytick',[]);
text(-0.1,1,'B','Units', 'Normalized','FontSize',12)

% dt = R.dt;
% step_tot = R.step_tot;
% seg_size = 4.5*10^3; % 2*10^4 for 2-pop, segmentation size for each plot
% seg_num = ceil(step_tot/seg_size);
% for i = 6
%     for seg = 6 %:seg_num
%         seg_ind = get_seg(step_tot, seg_size, seg);
%         t = (seg_ind-1)*dt; % ms
%         ax1 = subplot(2,2,2);
%         LFP_tmp = R.LFP.LFP_broad(i,seg_ind);
%         rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
%         plot(t, LFP_tmp,'b',t,rip_tmp,'k');
%         xlabel('Time(ms)','fontsize',10)
%         ylabel('LFP(a.u.)','fontsize',10)
%         legend('Broadband','Ripple-band')
%         text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
%         scales = R.LFP.wavelet.scales;
%         x_tmp = R.LFP.LFP_ripple(i,seg_ind);
%         coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
%         ax2 = subplot(2,2,4); % Scaleogram with pseudo-Frquency
%         freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
%         imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
%         xlabel('Time(ms)','fontsize',10)
%         ylim(freqrange)
%         ylabel('Frequency(Hz)','fontsize',10)
%         text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
%         % Link axes to synchronise them when zooming
%         linkaxes([ax1 ax2],'x');
%         set(ax1, 'xlim', minmax(t));
%     end
% end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig7 % this is the trick!!