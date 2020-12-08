function [V, tauy] = CollectTE(var,filenames)
% [result, loop_num] = CollectVectorYG(var, data)
%
% Collect data into a vector. The function automatically loops through any
% files with an extension "*RYG.mat" in the current working directory.
%
% var : the name of the variable that should be loaded from *RYG.mat files
% data: an expression to be evaluated to generated the data to be collected
%
% For example:
%   >> var = 'Analysis'
%   >> data = 'mean(Analysis.rate{1})'
%   >> [result, loop_num] = CollectCellYG(var, data)
%
% [result, loop_num] = CollectCellYG(var, data, filenames)
%
show_progress = 0;

% Prepare files
if nargin == 1
    dir_strut = dir('85*transferE.mat');
elseif nargin == 2
    dir_strut = dir(filenames);
end

num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

V = [];
tauy = [];
fprintf('Collecting data %s from %d files: \n',var, num_files);
for i = 1:num_files
    if  show_progress == 1
        fprintf('\t Loading data %s from file %s...',var, files{i});
    end
    load(files{i},var,'Opt');
    if  show_progress == 1
        fprintf('done.\n');
    end
    expr = sprintf('data_tmp = %s;', var);
    eval(expr);
    if isempty(data_tmp)
        warning('empty data')
    end
    V = [V; data_tmp ]; %#ok<AGROW>
    tauy = [tauy Opt.tauy];
       
    clear data_tmp; % clear it! Otherwise it could be misused by the consecutive loops.
    
end
[tauy,I] = sort(tauy,'descend');
V = V(I,:);
fprintf('\n');
end