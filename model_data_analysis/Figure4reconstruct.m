function Figure4reconstruct(varargin)
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
%
%   If the firing history is too long, data will be segmented into several
%   history fractions and plotted separately.
% % REFERENCE:Local generation and propagation of ripples along the septotemporal axis of the hippocampus
disp('plot_SWR...');
R = load('0001-201608271253-02473_in_1472266547975_out_RYG.mat');
R = get_grid_firing_centre(R);
load('0001-201608271255-46526_config_data.mat');

if nargin <= 1
    save_figure = 1; % default
end
if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end

dt = R.dt;
step_tot = R.step_tot;
seg_size = 10^3; % 2*10^4 for 2-pop, segmentation size for each plot
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_broad);
freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;
% Loop number for PBS array job
loop_num = 0;

% fig_n = 1;
for seg = 1:seg_num
    max_CData = [];
    min_CData = [];
    loop_num = loop_num + 1;
    
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    for i = 1:no
        seg_ind = get_seg(step_tot, seg_size, seg);
        t = (seg_ind-1)*dt*1e-3; % second
        x_tmp = R.LFP.LFP_ripple(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        min_C = min(transpose(coeffs_tmp));
        min_C = min(min_C(:));
        min_CData = [min_CData min_C];
        max_C = max(transpose(coeffs_tmp));
        max_C = max(max_C(:));
        max_CData = [max_CData max_C];
    end
    min_CData = min(min_CData(:));
    max_CData = max(max_CData(:));
    for i = 1:no
        seg_ind = get_seg(step_tot, seg_size, seg);
        t = (seg_ind-1)*dt*1e-3; % second
        x_tmp = R.LFP.LFP_ripple(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
        caxis([min_CData,max_CData]);
        colorbar('off');
        xlim([min(t) min(t) + seg_size/10^4]);
        ylim(freqrange);
        xlabel('s');
        ylabel('Hz');
        xlimits = get(gca,'xlim');
        set(gca,'Xtick',xlimits);
        ylimits = get(gca,'ylim');
        set(gca,'Ytick',ylimits);
        ax1 = gca;
        Position = get(ax1,'Position');
        ax2 = axes('Position',Position);
        set(ax2, 'visible', 'off')
        rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
        hold on;
        plot(t,rip_tmp, 'Parent',ax2,'Color','w','LineWidth',6);
        xlim(ax2,[min(t) min(t) + seg_size/10^4]);
        set(gcf,'color',[1 1 1]); % Set the figure frame color to white
        set(gca,'color',[1 1 1]); % Set the axis frame color to white
        name = sprintf(['%04d_SWR_%d_%d.pdf'],R.ExplVar.loop_num,i,seg);
        if save_figure == 1
            fprintf('\t Saving figure...');
            saveas(gca,name);
            fprintf('done.\n');
        else
            next = input('\t Next figure?');
            delete(gcf);
        end
    end
end
end
