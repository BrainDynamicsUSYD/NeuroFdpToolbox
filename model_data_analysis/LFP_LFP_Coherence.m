function LFP_LFP_Coherence(stdin)
% LFP-LFP coherence based on distance
if nargin == 0
    dir_strut = dir('*RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for i = 1:num_files
        files{i} = dir_strut(i).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end
[Lattice, ~] = lattice_nD(2, 1.5); % creat 4*4
dist = zeros(16);
window_size = 8; % ~2 Hz
for i = 1:16
    d = lattice_nD_find_dist(Lattice,1.5,i);
    dist(i,:) = d';
end
dist = tril(dist);
Dis = unique(sort(dist));
unit = 15.5*10; % um
for j = 2:length(Dis)
    [row,col] = find(dist == Dis(j));
    Cxy_mat = [];
    for i = 1:num_files % [1 12 17 18]
        % start form .mat files
        fprintf('Loading RYG.mat file %s...', files{i});
        R = load(files{i});
        disp('done.\n');
        x = R.LFP.LFP_gamma(row,10:10:end);
        x = x';
        y = R.LFP.LFP_gamma(col,10:10:end);
        y = y';
        [Cxy,F] = mscohere(x,y,[],[],[],1000);
        Cxy_mat = [Cxy_mat; nanmean(Cxy',1)];
        F = F';
    end
    n = round(length(F)/5); % Fmax=500Hz here
    subplot(5,1,j - 1)
    raw = nanmean(Cxy_mat(:,1:n),1);
    simple = tsmovavg(raw,'s',window_size,2);
    plot(F(1:n),simple)
    d = round(unit*Dis(j));
    ylabel(['d=',num2str(d),'um'])
end
xlabel({'Frequency(Hz)','LFP-LFP Cohernece'})
end