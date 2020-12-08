% test/watch mat by mat

dir_strut = dir('*RYG.mat');
% dir_strut2 = dir('*0_neurosamp.mat');
% dir_strut3 = dir('*1_neurosamp.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    reproduce_reference_no3F1A(R);
    power_spectrum(R);
    next = input('\t Next Mat?');
    close all;
end