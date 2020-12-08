function RasterPlotYL2(R,Neuron,sample_color,seg_ind)
% Sample color should be the size of sample neurons or empty
% Adjusting from function RasterPlotYL2.m
% Working for WM plotting, especially for overlapping items

% Input check and default values
seg = 1;
if nargin < 3
    sample_color = [];
    seg_ind = [];
end
seg_size = 4*10^4; % 4*10^4 for 2-pop, segmentation size for each plot 40s
sample_color_unit_range = 0;
% if 1, the default color range [0, 1] will be used
% if 0, sample_color vector will be scaled to [0, 1] automatically
text_fontsize = 10;
% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
% Segmetation
if isempty(seg_ind)
    seg_ind = get_seg(step_tot, seg_size, seg);
end
% Dump fields
num_spikes = R.reduced.num_spikes{1}(seg_ind);
spike_hist = R.reduced.spike_hist{1}(:,seg_ind);
T = (seg_ind-1)*dt;
Neuron = [Neuron{:}];
% Plot raster plot
if nnz(num_spikes) > 0
    if length(unique(Neuron)) == size(Neuron,1)*size(Neuron,2)
        ind_sample = Neuron(:); % sort(randperm(N(pop_ind),sample_size));
    else
        ind_sample = Neuron(:,1);
        ind_sample = [ind_sample;setdiff(Neuron(:,2),Neuron(:,3));intersect(Neuron(:,2),Neuron(:,3));setdiff(Neuron(:,3),Neuron(:,2))];
    end
    
    if isempty(sample_color)
        [Y,X,~] = find(spike_hist(ind_sample,:));
        xdata = ( [X(:)'; X(:)']+seg_ind(1)-1)*dt*1e-3; % sec
        ydata = [Y(:)'-1;Y(:)'];
        line(xdata, ydata,'Color','k');
    else
        hold on;
%         if sample_color_unit_range == 0 % automatically scale to [0,1]
%             sample_color = (sample_color-min(sample_color))/(diff(minmax(sample_color)));
%         end
        sample_color(sample_color == 0) = eps; % otherwise ceil() will give 0
        jetmap = colormap('parula(1000)');
        
        for i = 1:length(ind_sample)
            color_tmp = jetmap( ceil(sample_color(i)*1000),:);
            [Y,X,~] = find(spike_hist(ind_sample(i),:));
            xdata = ( [X(:)'; X(:)']+seg_ind(1)-1)*dt*1e-3; % sec
            ydata = [Y(:)'-1;Y(:)']+i-1;
            line(xdata, ydata, 'Color', color_tmp,'linewidth',0.1);
        end
        colormap('parula');
    end
    
    ylim([0,length(ind_sample)]);
    ylabel('Neurons','fontsize',text_fontsize)
    xlabel('Time(s)','fontsize',text_fontsize)
    xlim([T(1)-dt, T(1)+(length(seg_ind)-1)*dt]*1e-3); % make sure all the plots have the same axis scale
    
    % Keep tick lables while remove tick marks
%     set(gca,'ydir','reverse');
%     set(gca,'ytick',[],'ydir','reverse');
    box on;
end
end

