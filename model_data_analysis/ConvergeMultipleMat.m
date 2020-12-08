% read and converge multiple .mat to one
dir_strut = dir('*ratio0.05th20NTE.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
DC = zeros(1,30);
for i = 1:num_files
    R = load(files{i});
    DC(i) = R.dc;
end
dc = mean(vec2mat(DC,10),2)