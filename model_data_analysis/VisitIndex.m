function I = VisitIndex(Interval)% (varargin) % (R)
% two types spikes pattern: default ; modified
dir_strut = dir('*RYG.mat');
% dir_strut = dir('3DBurst3*minTime0SR1000.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
% % Loop number for PBS array job
% loop_num = 0;
I = zeros(1,num_files);

hw = 31;
% Interval = 30; % ms
Scale = 4; % 4*(7~10)um
pr = [];
[Lattice, ~] = lattice_nD(2, hw);

% coe = 2; % 2
% M = 20; % 20
for i = [1:13] % (num_files - 12):num_files
    %     loop_num = loop_num + 1;
    %     % For PBS array job
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    
%     R = get_grid_firing_centre(R,'mode','bayesian');
    t_mid = R.grid.t_mid;
    centre = R.grid.bayes.centre;
    factor = R.grid.bayes.bayes_factor_ln;
    centre(:,factor<log(100)) = NaN;
    
%     R = get_grid_firing_centreYL(R,'quick');
%     t_mid = R.grid.t_mid;
%     centre = R.grid.quick.centre;
%     width = R.grid.bayes.radius;
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
    
%     t_mid = cellfun(@(x) x(:,1)',R.WCentroids,'UniformOutput',false) ;
%     t_mid = t_mid{:};
%     centre = cellfun(@(x) x(:,2:3)',R.WCentroids,'UniformOutput',false) ;
%     centre = cell2mat(centre);
    
    for sites = 4:7
        for trials = 1:100
            Positions = 2*hw*rand(2,sites)-31;
            for t = 500:Interval:(length(t_mid)-Interval) % 500:Interval:(length(t_mid)-Interval)                
                candidates = centre(:,t:(t+Interval-1));
                
%                 candidates = NaN*ones(2,Interval);
%                 for k = 1:Interval
%                     if ismember(t+k-1,t_mid)
%                         candidates(:,k) = centre(:,t_mid==t+k-1);
%                     end
%                 end
                
                m = 0;
                for j = 1:sites
                    if min(Distance_xy(candidates(1,:),candidates(2,:),Positions(1,j),Positions(2,j),2*hw+1)) <= Scale
                        m = m + 1;
                    end
                end
                pr = [pr m/sites];
            end
        end
    end
    I(i) = sum(pr)/length(pr);
end
% index1 = sum(pr);
% index2 = sum(pr)/length(pr);
% disp(index2)
end