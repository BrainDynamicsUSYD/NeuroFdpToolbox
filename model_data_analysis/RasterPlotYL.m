function RasterPlotYL(R, pop_ind, seg,Neuron, sample_color,seg_ind)
% Sample color should be the size of sample neurons or empty
% Adjust from function raster_plot.m 
% No essential difference with default function, just to make it work in a
% way which I can handle well

% Input check and default values
if nargin < 3
    seg = 1;
end
if nargin < 5
    sample_color = [];
    seg_ind = [];
end

sample_size = 500; % sample neurons for raster plot
seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot

sample_color_unit_range = 0;
% if 1, the default color range [0, 1] will be used 
% if 0, sample_color vector will be scaled to [0, 1] automatically

text_fontsize = 12;

% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
N = R.N;
% dt = dt/1000; 

% Segmetation
if isempty(seg_ind)
    seg_ind = get_seg(step_tot, seg_size, seg);
end

% Dump fields
num_spikes = R.reduced.num_spikes{pop_ind}(seg_ind);
spike_hist = R.reduced.spike_hist{pop_ind}(:,seg_ind);
T = (seg_ind-1)*dt;

% Plot raster plot
if nnz(num_spikes) > 0
    % down-sampling
    if N(pop_ind) > sample_size
        ind_sample = Neuron; % sort(randperm(N(pop_ind),sample_size));
    else
        ind_sample = 1:1:N(pop_ind);
    end
    
    if isempty(sample_color)
        [Y,X,~] = find(spike_hist(ind_sample,:));
        xdata = ( [X(:)'; X(:)']+seg_ind(1)-1)*dt; % sec
        ydata = [Y(:)'-1;Y(:)'];
        line(xdata, ydata,'Color','k');
    else
        hold on;
        if sample_color_unit_range == 0 % automatically scale to [0,1]
            sample_color = (sample_color-min(sample_color))/(diff(minmax(sample_color)));
        end
        sample_color(sample_color == 0) = eps; % otherwise ceil() will give 0
        jetmap = colormap('parula(1000)');
        
        for i = 1:length(ind_sample)
            color_tmp = jetmap( ceil(sample_color(i)*1000),:);
            [Y,X,~] = find(spike_hist(ind_sample(i),:));
            xdata = ( [X(:)'; X(:)']+seg_ind(1)-1)*dt; % sec
            ydata = [Y(:)'-1;Y(:)']+i-1;
            line(xdata, ydata, 'Color', color_tmp,'linewidth',0.1);
        end
        colormap('parula');
    end    
    
    ylim([0,length(ind_sample)]);
    ylabel('Neurons','fontsize',text_fontsize)
    xlabel('Time(ms)')
    xlim([T(1)-dt, T(1)+(length(seg_ind)-1)*dt]); % make sure all the plots have the same axis scale
    
    % Keep tick lables while remove tick marks
    set(gca,'ytick',[],'ydir','reverse');
    box on;        
end
end

