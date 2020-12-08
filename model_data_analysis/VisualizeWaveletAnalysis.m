function VisualizeWaveletAnalysis
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    dt = R.dt;
    step_tot = R.step_tot;
    seg_size = 1e4; % 2s
    seg_num = ceil(step_tot/seg_size);
    [no, ~] = size(R.LFP.LFP_gamma);
    nos = 1:no;
    for i = nos
        for seg = 1:seg_num
            seg_ind = get_seg(step_tot, seg_size, seg);
            t = (seg_ind-1)*dt*1e-3; % second
            fs = 1/(dt*1e-3); % sampling frequency (Hz)
        end
    end
end
end