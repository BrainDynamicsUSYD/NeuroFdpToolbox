function D = DynamicsPatternDuration
% collect the spikes dynamics pattern duration
dir_strut = dir('Pattern0*.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
D = [];
for i = 7:13:num_files
    fprintf('Loading pattern ts %s...\n', files{i});
    R = load(files{i},'ts');
    Dif = round(diff(R.ts));
    Dif(Dif ~= 1) = 0;
    dur = diff([0 find(diff(Dif)) numel(Dif)]);
    if Dif(1) == 1
        start = 1;
    else
        start = 0;
    end
    duration = dur(start:2:end); % ms
    D = [D duration]; % ms
end
histogram(D)
xlabel('Duration(ms)')
ylabel('Count')
end