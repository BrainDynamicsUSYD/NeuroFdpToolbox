function distribution_subplot_1variable(stdin)
% A subplot format for distribution (1 variable)
% [v,loop_number] = CollectVectorYG('neuron_stats','neuron_stats.IE_ratio{1}');
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
figure(1)
for i = 1:6
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    subplot(2,3,i)
    histogram(R.neuron_stats.IE_ratio{1})
    xlabel(['decay_{GABA}=',num2str(2 + i)])
end
ax1 = axes('Position',[0 0 1 1],'Visible','off');% set nonvisible outside coordinate
axes(ax1);
title1 = 'IE ratio distribution';
text(0.34,0.97,title1,'fontsize',14)
saveas(gcf,'IE_ratio_distribution.pdf');
end