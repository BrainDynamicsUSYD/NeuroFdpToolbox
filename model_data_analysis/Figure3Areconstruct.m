function Figure3Areconstruct(varargin)
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
% REFERENCE:Local generation and propagation of ripples along the septotemporal axis of the hippocampus

disp('plot_SWR...');
R = load('0001-201608271253-02473_in_1472266547975_out_RYG.mat');

save_figure = 0;

if nargin < 1
    save_figure = 1; % default
end
if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end

dt = R.dt;
seg_size = 5*10^3; % 2*10^4 for 2-pop, segmentation size for each plot
freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;

i = 2;
start_ind = 151150;
seg_ind = start_ind:(start_ind + seg_size);
t = (seg_ind-1)*dt*1e-3; % second
x_tmp = R.LFP.LFP_ripple(i,seg_ind);
coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
colorbar;
xlim([min(t) min(t) + seg_size/10^4]);
ylim(freqrange);
xlabel('ms');
ylabel('Hz');
xlimits = get(gca,'xlim');
set(gca,'Xtick',xlimits);
set(gca,'XtickLabel',[0 500]);
ylimits = get(gca,'ylim');
set(gca,'Ytick',ylimits);
name = sprintf(['%04d_SWR_%d_%d.pdf'],R.ExplVar.loop_num,i,start_ind);
if save_figure == 1
    fprintf('\t Saving figure...');
    saveas(gca,name);
end
fprintf('done.\n');
end
