% IE ratio from sample mat
dir_strut0 = dir('*0_neurosamp.mat'); % E
dir_strut1 = dir('*1_neurosamp.mat'); % I
num_files = length(dir_strut0);
files0 = cell(1,num_files);
files1 = cell(1,num_files);
for i = 1:num_files
    files0{i} = dir_strut0(i).name;
    files1{i} = dir_strut1(i).name;
end
for i = 1:num_files % [1 12 17 18]
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files0{i});
    E = load(files0{i});
    I = load(files1{i});
    disp('done.\n');
    a = nanmean(E.I_GABA,2)./(nanmean(E.I_AMPA,2) + nanmean(E.I_ext,2));
    IE_ratio(i) = nanmean(a);
end