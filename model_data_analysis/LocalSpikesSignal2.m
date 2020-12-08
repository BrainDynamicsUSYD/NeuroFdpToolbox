function LocalSpikesSignal2(varargin)
% get local firing rate signal as "LFP"

% config = dir('0001*_config_data.mat');
% load(config.name,'LFP_centre_x','LFP_centre_y');
hw = 31;
LFP_centre_x = linspace(-hw, hw, 21); % E16:9  E100:21 E400:41
LFP_centre_y = linspace(-hw, hw, 21);
LFP_centre_x = LFP_centre_x(2:2:20); % E16(2:2:8)  E100(2:2:20) E400(2:2:40)
LFP_centre_y = LFP_centre_y(2:2:20);
[LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
LFP_centre_x = LFP_centre_x(:);
LFP_centre_y = LFP_centre_y(:);
dir_strut = dir('*_out_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
[Lattice, ~] = lattice_nD(2, hw);
S1 = 3.1; % 20 for cut off
Local1 = {};
for i = 1:length(LFP_centre_x)
    Local1{i} = find(lattice_nD_find_dist(Lattice,hw,LFP_centre_x(i),LFP_centre_y(i)) <= S1);
    
%     x1 = LFP_centre_x(i) - S1;
%     x2 = LFP_centre_x(i) + S1;
%     y1 = LFP_centre_y(i) - S1;
%     y2 = LFP_centre_y(i) + S1;
%     Local1{i} = find(Lattice(:,1)>=x1 & Lattice(:,1)<=x2 & Lattice(:,2)>=y1 & Lattice(:,2)<=y2);
end
% Loop number for PBS array job
loop_num = 0;
% k = [10 11];
for i = 1:num_files
    loop_num = loop_num + 1;   
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    for j = 1:length(LFP_centre_x)
        Spikes(j,:) = sum(R.spike_hist{1}(Local1{j}(:),:));
    end
    Spikes = reshape(full(Spikes),[100 50 R.step_tot/50]); % spikes in ms
    Spikes = squeeze(sum(Spikes,2));
%     Spikes = movsum(Spikes,[2 2],2); % window = 5 ms, gap = 1 ms
%     Spikes = Spikes(:,3:end-2);
    SpikesGrid = reshape(Spikes,[10 10 R.step_tot/50]);
    save([sprintf('%04g-',i),'LocalSpikesGridW5G5.mat'],'SpikesGrid')
end
end