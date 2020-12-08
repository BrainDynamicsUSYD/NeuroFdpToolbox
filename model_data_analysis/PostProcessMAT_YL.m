function PostProcessMAT_YL(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for i = 1:num_files
        files{i} = dir_strut(i).name;
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
for i = 1:num_files
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i}); % for performance, only load the necessary stuff
    
    disp('done.\n');
%     Figure4reconstruct(R);
    TimeFrequencyCImage(R);
%     R = get_grid_firing_centre(R);
%     grid = R.grid;
%     save(files{i},'grid', '-append');

%     sca = scatter(R);
%     s = sprintf('loop:%04i scatter:%d',R.ExplVar.loop_num,sca);
%     disp(s);
%     heatmap_single_bump_inR(R);% heatmap of the bump based on radius 
%     imagesc(flipud(vec2mat(R.Analysis.rate{1},63,63)));
%     colorbar
%     time = datestr(now,'yyyymmddTHHMMSS');
%     tit = sprintf('firing rate image of excitatory neurons,loop number = %04i time:%s',R.ExplVar.loop_num,time);
%     name = sprintf('%04i_image_firing_rate_%s.pdf',R.ExplVar.loop_num,time);
%     title(tit)
%     set(gcf,'renderer','zbuffer');
%     saveas(gcf,name);
    %%%%%%% do something here
    
    
    
    %     R_temp = rmfield(R_temp,'cluster');
    %     R_temp = cluster_sorted_rate(R_temp);
    %     cluster = R_temp.cluster;
    %     save(files{i},'cluster', '-append');
    
    
    % R_temp = get_CC_network(R_temp);
    % R_temp = get_CC_pop(R_temp, 1);
    % R_temp = get_CV2_ISI(R_temp);
    %     R_temp = get_ISI_low_high(R_temp);
    %     Analysis = R_temp.Analysis;
    %     save(files{i},'Analysis', '-append');
    
    
%     R_temp = avalanche_detect(R_temp);
%     avalanche = R_temp.avalanche;
%     save(files{i},'avalanche', '-append');
    
    
%     R_temp = get_grid_firing_centre(R_temp);
%     grid = R_temp.grid;
%     save(files{i},'grid', '-append');
    
%     R_temp = get_stPR(R_temp);
%     stPR = R_temp.stPR;
%     save(files{i},'stPR', '-append');
    

    %%R_temp = get_SWR(R_temp);
    %%SWR = R_temp.SWR;
    %%save(files{i},'SWR', '-append');
    %%plot_SWR({R_temp}, save_fig);
    
    
    % R_temp = rmfield(R_temp,{'C_rate','C_label','up_down_analysis'});
    
    
    %     R_temp = get_CC_pop(R_temp);
    %     R_temp = get_EI_current_crosscorr(R_temp);
    %     Balance = R_temp.Balance;
    %     Analysis = R_temp.Analysis;
    %     save(files{i},'Balance', 'Analysis', '-append');
    
    %     R_temp = cluster_sorted_rate(R_temp);
    %     cluster = R_temp.cluster;
    %     save(files{i},'cluster', '-append');
    %     ClusterRasterYG({R_temp}, save_fig);
    %     RasterYG({R_temp}, save_fig);
    
    %save(files{i},'-struct', 'R_temp', '-v7.3'); % -v7.3 for >2GB
    
    %     R_temp = get_neuron_sample_stats(R_temp);
    %     neuron_sample_stats = R_temp.neuron_sample_stats;
    %     save(files{i},'neuron_sample_stats', '-append');
    
    
    
    %%%%%%% do something above
    disp('Data processed and saved.');
    
end


end
