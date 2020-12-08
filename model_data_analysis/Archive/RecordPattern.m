function RecordPattern(varargin)
% adapt and base on visualize_grid_firing_centre.m function
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
num = [];
hw = 31;
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
        if ~isnan(x_centre(t)) && max(abs(R.grid.quick.centre(:,t))) < 31.5
            x_tmp = x_centre(t);
            y_tmp = y_centre(t);
            [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24))); % because win_len = 50 in get_grid_firing_centre.m
            spikingn = unique(spikingn);
            all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t)); % 81~1257
            incircle = sum(ismember(spikingn,all));
            %             if incircle/length(all) >= 0.07 && incircle > 23 % old setting: >=0.07 || >23
            if incircle/length(spikingn) > pi*width(t)^2/(2*hw+1)^2 && ( incircle/length(all) > 0.05 || incircle > 20 )
                x = [x x_tmp];
                y = [y y_tmp];
                ts = [ts t_mid(t)*0.1]; % ms
                r = [r width(t)];
                num = [num incircle];
            end
        end
    end
    for j = 1:length(x)-1
        d(j) = Distance_xy(x(j),y(j),x(j+1),y(j+1),2*hw+1);
        t(j) = ts(j+1) - ts(j);
    end
    v = d./t; % grid_d/ms:v*10um/ms mm/s
    save(['UPattern-',sprintf('%04g',loop_num),'.mat'],'x','y','r','ts','v','num');
end
end