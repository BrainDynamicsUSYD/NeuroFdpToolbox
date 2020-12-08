function cor = DurationVSFR % (R)
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
cor = zeros(1,num_files);
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    LFP_gamma = R.LFP.LFP{1}; % R.LFP.LFP_gamma; R.LFP.LFP{1}
    [no,~] = size(R.LFP.LFP_gamma); % R.LFP.LFP_broad; R.pop_stats.V_mean{1}
    ICI = [];
    Spikes = [];
    for i = 1:no
        [~,l1] = findpeaks(LFP_gamma(i,:));
        tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
        tempSpikes = zeros(1,length(l1)-1);
        for j = 2:length(l1)-1
            tempSpikes(j-1) = sum(R.num_spikes{1}(l1(j):l1(j+1)));
        end
        ICI = [ICI tempICI];
        Spikes = [Spikes tempSpikes];
    end
    % ind = find(ICI >= 5 & ICI <= 40);
    % ICI = ICI(ind);
    % Spikes = Spikes(ind);
    FR = Spikes./(ICI*1e-3)/R.N(1);
    ind = find(FR < 100);
    ICI = ICI(ind);
    FR = FR(ind);
    % dat = [ICI',FR'];
    % n = hist3(dat,[50 50]);
    % n1 = n';
    % n1(size(n,1) + 1, size(n,2) + 1) = 0;
    % xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    % yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
    % h = pcolor(xb,yb,n1);
    % set(h, 'EdgeColor', 'none');
    % h.ZData = ones(size(n1)) * -max(max(n));
    % colormap(hot)
    % oldcmap = colormap(gray);
    % colormap( flipud(oldcmap) );
    % colorbar
    [r,p] = corrcoef(ICI',FR');
    % if p(2) > 10^(-4)
    %     pvalue = ['p = ',num2str(p(2),'%.2f')];
    % else
    %     pvalue = ['p < 10^{-4}'];
    % end
    % %     title(['r = ',num2str(r(2),'%.2f'),pvalue])
    % tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
    % text(0.83,0.08,tp,'Units', 'Normalized','FontSize',8)
    % text(-0.28,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    % ylabel('Firing Rate(Hz)')
    % xlabel('Duration(ms)')
    cor(id_out) = r(2);
end
end