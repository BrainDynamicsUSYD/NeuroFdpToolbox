function d = PatternDuration
dir_strut = dir('Pattern*.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
for i = 7 % 1:num_files
    fprintf('Loading Pattern.mat file %s...\n', files{i});
    P = load(files{i});
    num = length(P.ts);
    binaryind = zeros(1,num-101);
    for j = 102:num-1
        d = Distance_xy(P.x(j),P.y(j),P.x(j-1),P.y(j-1),63);
        if P.ts(j)-P.ts(j-1)==1 && d < 18
            binaryind(j-100) = 1;
        end
    end
    ind = find([-1 diff(binaryind) -1] == -1);
    d = ind(2:end) - ind(1:end-1);
    histogram(d)
end
end