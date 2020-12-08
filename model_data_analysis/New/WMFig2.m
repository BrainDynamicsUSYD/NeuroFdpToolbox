% function WMFig2
figure_width = 7.8; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'WMFig2', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!
%%
figure
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
R = load(files{1});
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
LocalNeu = cell(1,R.ExplVar.NumP);
for i = 1:R.ExplVar.NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LocalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
LoalNeu = LocalNeu;
% LoalNeu = cell(1);
% LoalNeu{1} = LocalNeu{2};
for id = 32:num_files % normal; E; I % 11
subplot(7,1,[1,2])
R = load(files{id});
t_ind = 2.76e4:5e4;
try
    dt = R.reduced.dt;
    RasterPlotYL2(R,LoalNeu)
catch
    R.reduced.dt = 0.1;
    R.reduced.step_tot = length(R.spike_hist{1});
    R.reduced.spike_hist{1} = R.spike_hist{1};
    R.reduced.num_spikes{1} = full(sum(R.spike_hist{1},1));
    seg_ind = t_ind;
    RasterPlotYL2(R,LoalNeu,[],seg_ind)
end
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)

for i = 1:3
Signal = R.LFP{1}(i,:); % 1
fs = 1e4;
[wt,f,coi] = cwt(Signal,fs,'VoicesPerOctave',30);
for j = 1:length(coi)
    ind = find(f<=coi(j));
    wt(ind,j) = NaN;
end
ind = find(f<30 | f>140);
wt(ind,:) = [];
f(ind) = [];
tms = (0:numel(Signal)-1)/fs; % FRcertain
subplot(7,1,i+2)
surface(tms(t_ind),f,abs(wt(:,t_ind))) % cfs
colorbar off
caxis([0 12]) % [0 12]
axis tight
shading flat
xlabel('Time(s)','FontSize',10)
ylabel('Frequency(Hz)','FontSize',10)

% % Butterworth filter
% order = 4; % 4th order
% lowFreq = 30; % gamma band
% hiFreq = 140;
% Wn = [lowFreq hiFreq]/(fs/2);
% [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
% LFP_gamma = filter(b,a,Signal);
% lowFreq = 4; % gamma band
% hiFreq = 7;
% Wn = [lowFreq hiFreq]/(fs/2);
% [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
% LFP_theta = filter(b,a,Signal);
% LFP_gamma = (LFP_gamma-min(LFP_gamma))/(max(LFP_gamma)-min(LFP_gamma))*110+30;
% LFP_theta = (LFP_theta-min(LFP_theta))/(max(LFP_theta)-min(LFP_theta))*110+30;
% hold on
% plot3(tms(t_ind),LFP_gamma(t_ind),13*ones(1,length(t_ind)),'w') %,'LineWidth',1.5)
% hold on
% plot3(tms(t_ind),LFP_theta(t_ind),13*ones(1,length(t_ind)),'w') %,'LineWidth',1.5)
if i == 1
    text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
end
end
next = input('\t Next figure?');
close all
end
%%    
set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc WMFig2 % this is the trick!!
% end