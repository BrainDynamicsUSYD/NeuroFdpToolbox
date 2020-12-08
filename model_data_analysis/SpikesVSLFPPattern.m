% function SpikesVSLFPPattern
% visualize spikes pattern + LFP amp after Hilbert + space-temporal
% continuous burst

load('0015-201806291103-38833_in_1530234329865_0_neurosamp.mat','LFP_grid');
R = load('0015-201806291103-38833_in_1530234329865_out_RYG.mat');
load('3DBurst30015minTime30SR1000.mat','WCentroids')
RR.LFP.LFP{1} = reshape(LFP_grid,63^2,[]);
[RR] = GetBurst2(RR);
sigBinary = reshape(RR.LFP.GammaBurstEvent.is_burst,63,63,[]);
clear LFP_grid
LFP_grid = reshape(RR.LFP.LFP_gamma_hilbert_abs,63,63,[]);
clear RR
%%
showt = [];
for i = 1:13
    showt = [showt WCentroids{i}(:,1)'];
end
clim = minmax(reshape(LFP_grid,1,[]));
% Y = prctile(gamma_amp_grid(:),95);
S = full(R.spike_hist{1});
S = reshape(S,63,63,[]);
% vidObj = VideoWriter('0140SpikesLFP0msNormalization.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 30; % number of frames to display per second
% open(vidObj);
for i = 4291:10:4781 % (length(S)-1e4)
    LFP = LFP_grid(:,:,ceil(i/1));
    Spikes = sum(S(:,:,(i-25:i+24)+2e4),3);
    [row,col] = find(Spikes);
    imagesc(LFP',clim);
    hold on    
    [Row,Col] = find(sigBinary(:,:,i));
    if ismember(0.1*(i-1),showt)
         plot(Row,Col,'rs')         
    end
    hold on
    plot(row,col,'k.')
    ts = sprintf('time = %8.1f ms', i*0.1);
    title(ts);
    pause(0.1) % drawnow % pause(0.05)
    hold off
%     writeVideo(vidObj, getframe(gca));
end
% close(gcf);
% close(vidObj);
%%
S = full(R.spike_hist{1});
S = reshape(S,63,63,[]);
LFP_grid = reshape(R.LFP.LFP_gamma_hilbert_abs,63,63,[]);
% clear R
vidObj = VideoWriter('0013Spikes5msGammaLFPHilbertAbs.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20; % number of frames to display per second
open(vidObj);
for i = 10001:10:20001 % length(S) % 0.1 ms
    LFP = LFP_grid(:,:,ceil(i/10));
    Spikes = sum(S(:,:,(i-25:i+24)),3);
    [row,col] = find(Spikes);
    imagesc(LFP',clim);
    hold on
    plot(row,col,'k.')
    ts = sprintf('time = %8.1f ms', i*0.1);
    title(ts);
    pause(0.1) % drawnow % pause(0.05)
    hold off
    writeVideo(vidObj, getframe(gca));
end
close(gcf);
close(vidObj);
%% Pattern propagation snapshots baseline+2SD
fig = figure;
clim = minmax(reshape(LFP_grid(:,:,9e3:1e4),1,[]));
for i = 1:4
    t = [9080 9410 9700 9910];
    subplot(1,4,i)    
    LFP = LFP_grid(:,:,t(i));
    Spikes = sum(S(:,:,(t(i)-25:t(i)+24)+2e4),3);
    [row,col] = find(Spikes);
    imagesc(LFP',clim);
    hold on    
    [Row,Col] = find(sigBinary(:,:,t(i)));
%     if ismember(0.1*(t-1),showt)
%          plot(Row,Col,'rs','MarkerSize',2) %,'color',[1 0.8 0.8])         
%     end
%     hold on
    plot(row,col,'k.')
%     ts = sprintf('t = %8.1f ms',0.1*t);
%     title(ts);
    if i == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    end
end
%% Yifan's visualization
% function show_LFP_continous_and_grid(R)
dt = R.dt;

% down-sampling (skipped steps between video frames) for  video
down_sample = 10;

make_video = 1;
if make_video == 1
    writerObj = VideoWriter('LFPandSpikes2.avi');
    open(writerObj)
end

%%%%%% Prepare X (LFP-related data)
stamp = R.stamp;
samp_file = [stamp(1:end-3) '0_neurosamp'];
load(samp_file, 'gamma_power_grid'); % needs to be prepared in neurosamp.mat
load('3DBurstLFP0015minTime30SR1000P95.mat','WCentroids')
showt = [WCentroids{1}(:,1)' WCentroids{2}(:,1)' WCentroids{3}(:,1)'];
gamma_amp_grid = gamma_power_grid(:,:,1:10:end);
Y = prctile(gamma_amp_grid(:),95);
X = gamma_power_grid; % 3D matrix, 63x63xsteps
clear gamma_power_grid;
t_X = find(R.neuron_sample.t_ind{1}); % steps where gamma_power_grid is sampled

% check grid size, usually 63-by-63
[N_s, ~, steps] = size(X);
N = N_s^2;
if N ~= R.N(1)
    error('Not all of the excitatory neurons are sampled!')
end
fw = sqrt(N);
hw = (fw-1)/2;



%%%%%%% Prepare grid-lated data
ind_a_vec = R.grid.ind_ab(1,:);
ind_b_vec = R.grid.ind_ab(2,:);
[Lattice, ~] = lattice_nD(2, hw);
x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
t_grid = R.grid.t_mid; %   steps where get_grid_firing_centre is calculated
ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];
x_centre = R.grid.bayes.centre(1,:);
y_centre = R.grid.bayes.centre(2,:);
width = R.grid.bayes.radius;
is_pattern = R.grid.bayes.bayes_factor_ln > log(100);

%%%%%%% Prepare figure
figure('NumberTitle','off','Name','gamma Hilbert Power','color', 'w');
if make_video == 1
    set(gcf,'visible','off')
end
mm = minmax(reshape(X(:),1,[]));
xlim([-hw hw]);
ylim([-hw hw]);
hold on
h1 = imagesc(-hw:hw, -hw:hw,X(:,:,1), mm);
h2 = plot(0, 0, 'k.');
h3 = plot(0, 0, 'k.');
h4 = plot(0, 0, 'rs');
title('gamma Hilbert Power and fitted spiking pattern')

% start animation
t_start = 27; % 27
for i = t_start:down_sample: 3410 % steps
%     i/steps
    
    xlabel([num2str(t_X(i)*dt ),' ms']);
    
    
    pause(0.05)
    delete(h3)
    
    % Plot LFP-related data
    set(h1,'CData',X(:,:,i)');
    
    % Find the same step for grid and LFP data
    [t_diff,t_same] = min(abs(t_grid - t_X(i)));
    if t_diff ~= 0
        warning('Need to modify t_start to find match!')
    end
    
    % Plot spikes
    ind_range_tmp = ind_a_vec(t_same):ind_b_vec(t_same);
    set(h2,'XData',x_pos_o(ind_range_tmp),'YData',y_pos_o(ind_range_tmp));
    % Add spatial-temporal burst
    A = gamma_amp_grid(:,:,t_same-1999);
    GGrid = zeros(63);
    GGrid(A >= Y) = 1;
    [Row,Col] = find(GGrid);
    Col = Col - 32;
    Row = 64 - Row - 32;
    if ismember(t_same-1999,showt)
        set(h4,'XData',Row,'YData',Col);
    else
        set(h4,'XData',33,'YData',33);
    end
    
    % Plot detected spiking pattern if detected by Bayesian method
    if ~isnan(x_centre(t_same)) &&  is_pattern(t_same)
        x_tmp = x_centre(t_same);
        y_tmp = y_centre(t_same);
        r_cos = x_tmp+width(t_same)*cos(ang);
        r_sin = y_tmp+width(t_same)*sin(ang);
        h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
            r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
            r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
            r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
            r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
            r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
            r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
            r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
            r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
    else
        h3 = plot(0, 0, 'k.');
    end
    
    if make_video == 1
        writeVideo(writerObj, getframe(gcf));
    end
    
end
if make_video == 1
    close(writerObj);
end

% end
% end