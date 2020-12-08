function RecordPattern0(varargin)
% adapt from function RecordPattern.m function
% Collect data on default setting for pattern capture and Analysis
% location (x,y)
% corresponding time t
close all;
clc;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

%%% Loop number for PBS array job %%%
loop_num = 0;

x = [];
y = [];
ts = [];
r = [];
n = [];
hw = 31;
interval = 1; % ms; because the win_gap = 10 steps
[Lattice, ~] = lattice_nD(2, hw);
for i = 1:num_files
    
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    R = get_grid_firing_centre(R);
    t_mid = R.grid.t_mid;
    x_centre = R.grid.quick.centre(1,:);
    y_centre = R.grid.quick.centre(2,:);
    width = R.grid.quick.radius;
    
    for t = 1:length(t_mid)
        if ~isnan(x_centre(t)) % && max(abs(R.grid.quick.centre(:,t))) < 31.5
            x_tmp = x_centre(t);
            y_tmp = y_centre(t);
            [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24))); % because win_len = 50 in get_grid_firing_centre.m
            %             spikingn = unique(spikingn);
            %             all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t)); % 81~1257
            %             incircle = sum(ismember(spikingn,all));
            %             if incircle/length(spikingn) > pi*width(t)^2/(2*hw+1)^2 && ( incircle/length(all) > 0.1 || incircle > 20 )
            x = [x x_tmp];
            y = [y y_tmp];
            ts = [ts t_mid(t)*R.dt]; % ms
            r = [r width(t)];
            n = [n length(spikingn)];
            %                 num = [num incircle];
            %             end
        end
    end
    j = 1;
    count = 1;
    num = n(1);
    du = 1;
    while j < length(x)
        if ts(j+1)<ts(j)+interval+1 && Distance_xy(x(j),y(j),x(j+1),y(j+1),2*hw+1)<r(j)+r(j+1) % interval 1ms
            num = num + n(j+1);
            du = du + ts(j+1)-ts(j);
        else
            nums(count) = num;
            duration(count) = du;
            count = count + 1;
            num = n(j+1);
            du = 1;
        end
        j = j + 1;
    end
    save(['Pattern0-',sprintf('%04g-interval%d',loop_num,interval),'ms.mat'],'x','y','r','ts','nums','duration');
end
end