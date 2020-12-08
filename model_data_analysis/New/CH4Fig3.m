function CH4Fig3
figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'CH4Fig3', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% dir_strut2 = dir('*_config_data.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
id = 80; % 1;
R = load(files{id});
% R.reduced.step_tot = 7e3;
% R.reduced.num_spikes{1} = R.reduced.num_spikes{1}(1:)
% load(files2{1},'StiNeu')
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
LoalNeu = cell(1,R.ExplVar.NumP);
for i = 1:R.ExplVar.NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
subplot(3,1,3)
RasterPlotYL2(R,LoalNeu)
text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
for i = 1:2
    subplot(3,1,i)
    switch i
        case 1
            Signal = R.LFP.LFP{1}(1,:); % 1
            fs = 1e4;
            text(-0.1,1,'A','Units', 'Normalized','FontSize',12)
        case 2
            try
                FR = vec2mat(R.num_spikes{1},10);
            catch
                FR = vec2mat(full(sum(R.spike_hist{1})),10);
            end
            Signal = sum(FR,2)'/3969*1e3;
            fs = 1e3;
            text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
    end
%     Signal = Signal(fs:end);
    [wt,f,coi] = cwt(Signal,fs,'VoicesPerOctave',30);
    for j = 1:length(coi)
        ind = find(f<=coi(j));
        wt(ind,j) = NaN;
    end
    ind = find(f<30 | f>80);
    wt(ind,:) = [];
    f(ind) = [];
    tms = (0:numel(Signal)-1)/fs; % FRcertain
    surface(tms,f,abs(wt)) % cfs
%     set(gca,'YAxisLocation','right')
    colorbar off
    caxis([0 12]) % [0 12]
    axis tight
    shading flat
    xlabel('Time(s)','FontSize',10)
    ylabel('Frequency(Hz)','FontSize',10)
end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig3 % this is the trick!!
end