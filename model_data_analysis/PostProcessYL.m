function PostProcessYL(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory
% adjust from PostProcesYG.m function

% Prepare files
if nargin == 0
    dir_strut = dir('*out.h5');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end

% save figures
save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for id_out = 1:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = ReadH5_YL( files(id_out) ); % read .ygout file into matlab data struct
    
    %     meanx = R{1}.syn_sample.x_mean;
    %     flat = min(meanx(1:1.9e4));
    %     period = meanx(3e4:end);
    %     [pks,locs] = findpeaks(-period);
    %     locs = locs(-pks < (flat - 0.2));
    %     if isempty(locs) || sum((diff(locs) > 1e2) & (diff(locs) < 4e3)) < 1
    %         fprintf('Not good, terminating...\n')
    %         continue
    %     end
    
    R = AnalyseYG(R); % do some simple analysis
    fprintf('Starting getting Gamma...\n');
    R = GetGamma(R{1});
    R = get_grid_firing_centre(R);
%     R = get_grid_firing_centre(R,'mode','bayesian');
    R = GetBurst(R); % check it each time when using
    R = {R};
    SaveRYG(R);
    disp('Done');
    for ind = 1:length(R)
        if isfield(R{ind}, 'samp_file')
            Read_and_save_YGSamp(R{ind}.samp_file, R{ind});
        end
    end
    RasterYL(R, save_fig); % generate raster plot for spiking history
end
end

