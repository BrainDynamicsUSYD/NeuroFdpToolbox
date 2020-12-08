function [ OutData ] = ReadSampleH5( files )
% files are cell of array of char array (strings) defining the path of input files
% If given no argument, it searches for matches under CURRENT directory
% Only read I_AMPA and I_GABA in sample.h5
disp('ReadH5...');
tic;


% Prepare filenames
if nargin == 0
    % If given no argument, search for matches under CURRENT directory
    dir_strut = dir('*out.h5');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
else
    % all the corresponding input files
    files_dir = fileparts(files{1}); % assuming all the input files are under the same directory
    if ~isempty(files_dir) % if not under current directory
        in_all_dir = dir(strcat(files_dir,'/*out.h5'));
        in_all = cell(1,1);
        for id_in = 1:length(in_all_dir)
            in_all{id_in} = strcat(files_dir,'/',in_all_dir(id_in).name);
        end
    else % if under current directory
        in_all_dir = dir('/*out.h5');
        in_all = cell(1,1);
        for id_in = 1:length(in_all_dir)
            in_all{id_in} = in_all_dir(id_in).name;
        end
    end
end
OutData = cell(1,length(files));
if ~isempty(files)
    % Read ygout file(s) specified by "name" and write data to "OutData"
    for id_out = 1:length(files)
        % OutData{id_out}.file = files{id_out};
        [file_dir, file_name, ~] = fileparts(files{id_out});
        if ~isempty(file_dir)
            OutData{id_out}.stamp = strcat(file_dir,'/', file_name);
        else
            OutData{id_out}.stamp = file_name;
        end
        fprintf('Current ReadYG file is: %s\n', files{id_out});
        % neuron sample file
        stamps = file_name(1:end-3);
        for pop_ind = 1:2
            samp_file = [stamps num2str(pop_ind-1) '_neurosamp.h5'];
            samp_file_mat = [samp_file(1:end-2) 'mat'];
            if exist(samp_file,'file') == 2 && exist(samp_file_mat,'file') ~= 2 % 2 for .mat file
                fprintf('   Generating %s...', samp_file_mat);
                I_AMPA = try_h5read( samp_file,  '/I_AMPA' );
                I_GABA = try_h5read( samp_file,  '/I_GABA' );
                %                 I_K = try_h5read( samp_file,  '/I_K' );
                %                 I_ext = try_h5read( samp_file,  '/I_ext' );
                %                 V = try_h5read( samp_file,  '/V' );
                %                 I_leak = try_h5read( samp_file,  '/I_leak' );
                save(samp_file_mat, 'I_AMPA', 'I_GABA','-v7.3');
                clear 'I_AMPA'  'I_GABA' ;
                %                 save(samp_file_mat, 'I_AMPA', 'I_GABA','I_K','I_ext' ,'V', 'I_leak','-v7.3');
                %                 clear 'I_AMPA'  'I_GABA' 'I_K' 'I_ext' 'V''I_leak';
                fprintf('done\n');
            end
        end
    end
end
end