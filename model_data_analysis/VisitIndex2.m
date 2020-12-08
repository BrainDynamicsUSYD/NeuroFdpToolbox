function [Duration,Displacement] = VisitIndex2(varargin) % (R)
% two types spikes pattern: default ; modified
% recording time/displacement in a complete visiting of all sites
tic;
dir_strut = dir('3DBurst30*minTime0SR1000.mat'); % 3DBurst30*minTime0SR1000 % *RYG
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
% Loop number for PBS array job
loop_num = 0;
hw = 31;
Scale = 6; % 4*(7~10)um
% [Lattice, ~] = lattice_nD(2, hw);

% coe = 2; % 2
% M = 20; % 20

for sites = [10 12] % [4:7 10 12]
    Duration = cell(1,190); % num_files
    Displacement = cell(1,190);
    %     Interval = cell(1,130);
    loop_num = loop_num + 1;
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    for i = 131:num_files % (num_files - 12):num_files
        
        fprintf('Loading RYG.mat file %s...', files{i});
        R = load(files{i});
        disp('done.\n');
        
        % for spikes pattern
        %         t_mid = R.grid.t_mid;
        %         centre = R.grid.quick.centre;
        %         centre = R.grid.bayes.centre;
        %         factor = R.grid.bayes.bayes_factor_ln;
        %         centre(:,factor<log(100)) = NaN;
        
        % for LFP pattern
        WCentroids = cell2mat(R.WCentroids');
        t_mid = WCentroids(:,1);
        centre = WCentroids(:,2:3)';
        
        %     t_mid = R.grid.t_mid;
        %     centre = R.grid.quick.centre;
        %     width = R.grid.quick.radius;
        
        %     for t = 1:length(t_mid)
        %         if ~isnan(centre(1,t))
        %             x_tmp = centre(1,t);
        %             y_tmp = centre(2,t);
        %
        %             % adding modification algorithm as criteria for proper spike pattern
        %             [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        %             spikingn = unique(spikingn);
        %             all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        %             incircle = sum(ismember(spikingn,all));
        %             if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
        %             else
        %                 centre(:,t) = NaN;
        %             end
        %         end
        %     end
        
        NoOverlap = 2*ones(1,sites);
        for trials = 1:100
            while sum(NoOverlap) > sites
                Positions = 2*hw*rand(2,sites)-31;
                for j = 1:sites
                    NoOverlap(j) = sum(lattice_nD_find_dist(Positions',hw,j) < 2*Scale);
                end
            end
            IndM = zeros(sites,length(t_mid));
            startflag = 0;
            %             tstart2 = 0;
            periodind = zeros(1,sites);
            for j = 1:sites
                IndM(j,:) = (Distance_xy(centre(1,:),centre(2,:),Positions(1,j),Positions(2,j),2*hw+1) <= Scale);
            end
            for t = 1:length(t_mid)-1 % 500:Interval:(length(t_mid)-Interval)
                if startflag == 0 && sum(IndM(:,t)) > 0
                    startflag = 1;
                    periodind = zeros(1,sites);
                    displacement = 0;
                    tstart = t_mid(t); % t_mid(t); % t;
                    %                     if tstart2 ~= 0
                    %                         interval = tstart - tstart2;
                    %                         Interval{i} = [Interval{i} interval];
                    %                     end
                end
                if startflag == 1
                    ind = find(IndM(:,t));
                    periodind(ind) = 1;
                    dis = Distance_xy(centre(1,t),centre(2,t),centre(1,t+1),centre(2,t+1),2*hw+1);
                    if ~isnan(dis)
                        displacement = displacement + dis;
                    end
                end
                if sum(periodind) == sites
                    duration = t_mid(t) - tstart; % t_mid(t) - tstart; % t - tstart;
                    Duration{i-130} = [Duration{i-130} duration]; % i-130
                    Displacement{i-130} = [Displacement{i-130} displacement];
                    startflag = 0;
                    periodind = zeros(1,sites);
                    %                     tstart2 = t_mid(t);
                end
            end
        end
    end
    save(['0131-0320LFPVisitIndexTwoNoOverlap',sprintf('%d',sites),'SitesScale5.mat'],'Duration','Displacement');
    toc;
end
end