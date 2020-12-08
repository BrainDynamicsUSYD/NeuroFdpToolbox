% [rate,~] = CollectVectorYG('Analysis','Analysis.rate');
% [nn,loop] = CollectVectorYG('LFP','[LFP.ripple_event.ripple_start_steps{:}]');
function read_collect_R(stdin)
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

% mean_rate_E = zeros(1,9);
% mean_CVISI_E = zeros(1,9);
% % mean_rate_I = zeros(1,9);
% mean_CVISI_I = zeros(1,9);
% j = 1;
for i = [1:22]
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    mean(R.neuron_stats.IE_ratio{1})
%     mean_rate_E(i) = nanmean(R.Analysis.rate{1});
%     mean_CVISI_E(j) = nanmean(sqrt(R.Analysis.CV2_ISI{1}));
% %     mean_rate_I(i) = nanmean(R.Analysis.rate{2});
%     mean_CVISI_I(j) = nanmean(sqrt(R.Analysis.CV2_ISI{2}));
%     j = j + 1;
end
% figure(1)
% subplot(2,1,1)
% % imagesc(vec2mat(mean_rate_E,3));
% % colorbar
% decay = [3:6 8 10 14:3:20];
% bar(decay,mean_CVISI_E)
% % set(gca,'XTick',[1:9])
% % set(gca,'XTickLabel',[3:6 8 10 14:3:20])
% % set(gca,'YTick',[1:3])
% % set(gca,'YTickLabel',[0.9:0.1:1.1])
% xlabel('decay_{GABA}(ms)')
% % ylabel('Phi_E')
% % title('firing rate(E)')
% % figure(2)
% % imagesc(vec2mat(mean_CVISI_E,3));
% % colorbar
% % set(gca,'XTick',[1:3])
% % set(gca,'XTickLabel',[0.9:0.1:1.1])
% % set(gca,'YTick',[1:3])
% % set(gca,'YTickLabel',[0.9:0.1:1.1])
% % xlabel('Phi_I')
% % ylabel('Phi_E')
% title('CV_{ISI}(E)')
% % figure(3)
% % imagesc(vec2mat(mean_rate_I,3));
% % colorbar
% subplot(2,1,2)
% bar(decay,mean_CVISI_I)
% % set(gca,'XTick',[1:18])
% % set(gca,'XTickLabel',[3:6 8 10 14:3:20])
% % set(gca,'YTick',[1:3])
% % set(gca,'YTickLabel',[0.9:0.1:1.1])
% xlabel('decay_{GABA}(ms)')
% % ylabel('Phi_E')
% % title('firing rate(I)')
% % figure(4)
% % imagesc(vec2mat(mean_CVISI_I,3));
% % colorbar
% % set(gca,'XTick',[1:3])
% % set(gca,'XTickLabel',[0.9:0.1:1.1])
% % set(gca,'YTick',[1:3])
% % set(gca,'YTickLabel',[0.9:0.1:1.1])
% % xlabel('Phi_I')
% % ylabel('Phi_E')
% title('CV_{ISI}(I)')
% saveas(gcf,'CV_ISI.pdf');
end