function input_rhythm(stdin)
% input rhythm
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
for i = [1:9]
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    t = (1:1e3)*0.1;
    t1 = (5e4 + 1):(5e4 + 1e3);
    subplot(3,3,i)
    plot(t,R.pop_stats.I_input_mean{1}(t1))
end
end