function LFPCrossCovariance(stdin)
% intermittent LFP Cross-Covariance
tic;
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
% window_size = 8; % ~2 Hz
for i = 1:16
    d = lattice_nD_find_dist(Lattice,1.5,i);
    dist(i,:) = d';
end
dist = tril(dist);
% unit = 15.5*10; % um
[row,col] = find(dist == 1);
window = 500; % 50 ms
Cxy_mat = [];
for i = 1%:num_files % [1 12 17 18]
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    for j = (1e4 + 1):(1e4 + 3e4)
        c = [];
        for k = 1:length(row)
            x = R.LFP.LFP_gamma(row(k),j:(j+window));
            y = R.LFP.LFP_gamma(col(k),j:(j+window));
            [cc,lags] = xcov(x,y,100); % maxlag = 10 ms
            c = [c;cc];
        end
        Cxy_mat = [Cxy_mat nanmean(c,1)'];
    end
end
% raw = nanmean(Cxy_mat(:,1:n),1);
% simple = tsmovavg(raw,'s',window_size,2);
imagesc(Cxy_mat)
colorbar
toc;
% ylabel(['d=',num2str(d),'um'])
end