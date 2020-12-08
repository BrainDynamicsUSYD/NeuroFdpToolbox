function SpikesPattern(varargin)
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
A = [];
B = [];
C = [];
% % Loop number for PBS array job
% loop_num = 0;
for id_out = 1 % :10 % 5:13:130 % num_files
    %      % For PBS array job
    %     loop_num = loop_num + 1;
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    if ~isfield(R,'grid') || ~isfield(R.grid,'bayes')
        R = get_grid_firing_centre(R,'mode','bayesian');
        SaveRYG_nocell(R)
    end
    
    %     raw = ~isnan(R.grid.quick.radius);
    raw = ~isnan(R.grid.bayes.radius); % .grid
    ind = find(R.grid.bayes.bayes_factor_ln <= log(100)); % .grid
    raw(ind) = 0;
    
    % coe = 2;
    % M = 20;
    tempConti = 10; % ms
    hw = 31;
    t_mid = R.grid.t_mid;
    x_centre = R.grid.bayes.centre(1,:);
    y_centre = R.grid.bayes.centre(2,:);
    [Lattice, ~] = lattice_nD(2, hw);
    width = R.grid.bayes.radius;
    i = 1;
    raw2 = [];
    while i <= length(raw)-tempConti
        switch raw(i)
            case 1
                switch raw(i+1)
                    case 1
                        if Distance_xy(x_centre(i),y_centre(i),x_centre(i+1),y_centre(i+1),2*hw+1) > sum(width([i,i+1])) % OR:max
                            raw2 = [raw2 1 NaN];
                        else
                            raw2 = [raw2 1];
                        end
                        i = i + 1;
                    case 0
                        ind = find(raw(i+1:i+tempConti));
                        if ~isempty(ind)
                            if Distance_xy(x_centre(i),y_centre(i),x_centre(i+ind(1)),y_centre(i+ind(1)),2*hw+1) <= sum(width([i,i+ind(1)]))
                                raw2 = [raw2 ones(1,ind(1))];
                            else
                                raw2 = [raw2 raw(i:i+ind(1)-1)];
                            end
                            i = i + ind(1);
                        else
                            raw2 = [raw2 raw(i:i+tempConti)];
                            i = i + tempConti + 1;
                        end
                end
            case 0
                raw2 = [raw2 0];
                i = i + 1;
        end
    end
    
    MYM = zeros(size(raw));
    for t = 1:length(t_mid)
        %         if raw(t) == 1 % ~isnan(x_centre(t))
        %             x_tmp = x_centre(t);
        %             y_tmp = y_centre(t);
        % adding modification algorithm as criteria for proper spike pattern
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-5):(t_mid(t)+4)));
        %             spikingn = unique(spikingn);
        %             all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        %             incircle = sum(ismember(spikingn,all));
        MYM(t) = length(spikingn); % incircle % length(spikingn)
        %             if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
        %                 MYM(t) = 1;
        %             end
        %         end
    end
    %     raw = raw.*MYM;
    
    raw = raw2;
    raw(isnan(raw)) = 0;
    [~,high_du1,~,start,~] = seq_postprocess(raw,1); % ms
    raw = raw2;
    raw(isnan(raw)) = [];
    [~,~,low_du1,~,~] = seq_postprocess(raw,1); % ms
    S = zeros(1,length(high_du1)); % size
    indraw2 = start(1);
    indraw = start(1);
    i = 1;
    raw2(isnan(raw2)) = 2;
    while indraw2 <= length(raw2) && i <= length(S)
        switch raw2(indraw2)
            case 1
                S(i) = S(i) + MYM(indraw);
                indraw = indraw + 1;
            case 2 % represent NaN
                i = i + 1;
            case 0
                if raw2(indraw2-1) ~= 0
                    i = i + 1;
                end
                indraw = indraw + 1;
        end
        indraw2 = indraw2 + 1;
    end
    %     for i = 1:length(high_du1)
    %         S(i) = sum(MYM(start(i):(start(i)+high_du1(i)-1)));
    %     end
    C = [C S];
    A = [A high_du1];
    B = [B low_du1];
end
% save('BayesFactorAllSize.mat','A')
fprintf('Saving done...\n')
%% version2 for spikes pattern
for id_out = 1 % :10 % 5:13:130 % num_files
    %      % For PBS array job
    %     loop_num = loop_num + 1;
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    if ~isfield(R,'grid') || ~isfield(R.grid,'bayes')
        R = get_grid_firing_centre(R,'mode','bayesian');
        SaveRYG_nocell(R)
    end
    
    %     raw = ~isnan(R.grid.quick.radius);
    raw = ~isnan(R.grid.bayes.radius); % .grid
    ind = find(R.grid.bayes.bayes_factor_ln <= log(100)); % .grid
    raw(ind) = 0;
    
    % coe = 2;
    % M = 20;
    tempConti = 10; % ms
    hw = 31;
    t_mid = R.grid.t_mid;
    x_centre = R.grid.bayes.centre(1,:);
    y_centre = R.grid.bayes.centre(2,:);
    [Lattice, ~] = lattice_nD(2, hw);
    width = R.grid.bayes.radius;
    i = 1;
    count = 1;
    Spikespattern = cell(count);
    while i <= length(raw)-tempConti
        switch raw(i)
            case 1
%                 switch raw(i+1)
%                     case 1
                        Spikespattern{count} = [Spikespattern{count};t_mid(i) x_centre(i) y_centre(i) width(i)];
%                         if Distance_xy(x_centre(i),y_centre(i),x_centre(i+1),y_centre(i+1),2*hw+1) > sum(width([i,i+1])) % OR:max
%                             count = count + 1;
%                             Spikespattern{count} = [];
%                         end
                        i = i + 1;
%                     case 0
%                         ind = find(raw(i+1:i+tempConti));
%                         if ~isempty(ind)
%                             if Distance_xy(x_centre(i),y_centre(i),x_centre(i+ind(1)),y_centre(i+ind(1)),2*hw+1) <= sum(width([i,i+ind(1)]))
%                                 Spikespattern{count} = [Spikespattern{count};t_mid(i:i+ind(1)-1)' x_centre(i:i+ind(1)-1)' y_centre(i:i+ind(1)-1)' width(i:i+ind(1)-1)'];
%                             else
%                                 Spikespattern{count} = [Spikespattern{count};t_mid(i) x_centre(i) y_centre(i) width(i)];
%                                 count = count + 1;
%                                 Spikespattern{count} = [];
%                             end
%                             i = i + ind(1);
%                         else
%                             Spikespattern{count} = [Spikespattern{count};t_mid(i) x_centre(i) y_centre(i) width(i)];
%                             count = count + 1;
%                             Spikespattern{count} = [];
%                             i = i + tempConti + 1;
%                         end
%                 end
            case 0
                i = i + 1;
                if ~isempty(Spikespattern{count})
                    count = count + 1;
                    Spikespattern{count} = [];
                end
        end
    end
end
% save('BayesFactorAllSize.mat','A')
fprintf('Saving done...\n')
%% Spikes pattern propagation snapshots
fig = figure;
tsteps = 3.016e3;
ind_a_vec = R.grid.ind_ab(1,:);
ind_b_vec = R.grid.ind_ab(2,:);
hw = 31;
fw = 2*hw+1;
ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];
x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
for i = 1:8
    t = (tsteps-16)/10 + 5*i;
    subplot(2,4,i)
    axis equal;
    box on;
    set(gca,'xtick',[],'ytick',[]);    
    xlim([-hw hw]);
    ylim([-hw hw]);
    hold on;
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    x_tmp = x_centre(t);
    y_tmp = y_centre(t);
    r_cos = x_tmp+width(t)*cos(ang);
    r_sin = y_tmp+width(t)*sin(ang);
    hold on
    plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
        r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
        r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
        r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
        r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
        r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
        r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
        r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
        r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
    hold on
    plot(x_tmp,y_tmp, 'r>', 'MarkerSize', 8);
    ts = sprintf('t = %8.1f ms',t);
    title(ts);
    if i == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    end
end

%% plot distribution & fit levy alpha stable
f = fitdist([A(:);-A(:)],'stable')
f2 = fitdist([A(:)],'stable')
figure
[N,edge]=histcounts(A,51,'normalization','pdf');
subplot(2,2,1)
shadedErrorBar(edge(2:end),N,zeros(size(N)))
title('Size') % Interval % Size
xlabel('Size') % Duration(ms) % number of spikes
ylabel('Probability')
set(gca,'yscale','log','xscale','log')
set(gca,'linewidth',1.5,'fontsize',15)

subplot(2,2,2)
[N,edge]=histcounts(A(:),51,'normalization','pdf');
loglog(edge(2:end),N,'.')
title('Size')
xlabel('Size')
ylabel('Probability')
set(gca,'linewidth',1.5,'fontsize',15)

subplot(2,2,3)
m = max(A);
x = linspace(-m,m,1000);
plot(x,pdf(f,x),'r-','linewidth',2)
hold on
h=histogram([A(:);-A(:)],'normalization','pdf');% #bins:60
xlabel('Size')
ylabel('Probability')
set(gca,'linewidth',1.5,'fontsize',15)

subplot(2,2,4)
x = linspace(min(A),m,1000);
plot(x,pdf(f2,x),'r-','linewidth',2)
hold on
h=histogram([A(:)],'normalization','pdf');% #bins:60
xlabel('Size')
ylabel('Probability')
set(gca,'linewidth',1.5,'fontsize',15)
% saveas(gcf,'Window10ms.eps')

% %% add break time tolerance
% btime = 3;
% ind = find(low_du1 <= btime);
% Raw = raw;
% for i = ind
%     Raw(low_start1(i):(low_start1(i)+low_du1(i)-1)) = 1;
% end
% [~, high_du2, low_du2, ~, ~] = seq_postprocess(Raw,1); % ms
% %% add jump step size threshold
% Raw = raw;
% jumpss = 15;
% IND = find(R.grid.quick.jump_dist > jumpss);
% modify = [];
% for i = 1:length(Raw)
%     modify = [modify Raw(i)];
%     if ismember(i,IND)
%         modify = [modify 0];
%     end
% end
% [~, high_du3, low_du3, ~, ~] = seq_postprocess(modify,1); % ms
end

function SaveRYG_nocell( Result )
%This function saves Result_cell as explanatory and response
%variables to .mat file
%   Detailed explanation goes here

name_temp = strcat( Result.stamp, '_RYG.mat');
fprintf('Saving results into file %s...\n', name_temp);
save(name_temp, '-struct', 'Result', '-v7.3'); % -v7.3 for >2GB

disp('Saving done.')

end