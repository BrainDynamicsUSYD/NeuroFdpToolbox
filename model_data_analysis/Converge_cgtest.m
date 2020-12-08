function [X1,X2,Y1,Y2] = Converge_cgtest
% read and converge multiple .mat to one
% specifically working in directory change_gbalance
dir_strut = dir('teen*.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
X1 = [];
X2 = [];
Y1 = [];
Y2 = [];

for i = 66:num_files
    R = load(files{i});
    X1 = [X1 nanmean(R.nX1)];
    X2 = [X2 nanmean(R.nX2)];
    Y1 = [Y1 nanmean(R.nY1)];
    Y2 = [Y2 nanmean(R.nY2)];
end

X1 = vec2mat(X1,13);
X2 = vec2mat(X2,13);
Y1 = vec2mat(Y1,13);
Y2 = vec2mat(Y2,13);

X1 = nanmean(X1);
X2 = nanmean(X2);
Y1 = nanmean(Y1);
Y2 = nanmean(Y2);
end