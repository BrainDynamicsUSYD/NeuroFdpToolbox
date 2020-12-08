function WMFig1
figure_width = 15.6; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'WMFig1', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!
%%
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
R = load(files{1});
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
LoalNeu = cell(1,R.ExplVar.NumP);
for i = 1:R.ExplVar.NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
% text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
id = [11+21*93,1+21*11,21]; % normal; E; I
NumP = R.ExplVar.NumP;
bin = 500;
for i = 1:3
%     subplot(3,3,3*i-2)
    R = load(files{id(i)});
%     try
%         dt = R.reduced.dt;
%         RasterPlotYL2(R,LoalNeu)
%     catch
%         R.reduced.dt = 0.1;
%         R.reduced.step_tot = length(R.spike_hist{1})-2e4;
%         R.reduced.spike_hist{1} = R.spike_hist{1}(:,1:end-2e4);
%         R.reduced.num_spikes{1} = full(sum(R.spike_hist{1}(:,1:end-2e4),1));
%         seg_ind = get_seg(R.reduced.step_tot, 4e5, 1);
%         RasterPlotYL2(R,LoalNeu,[0.1*ones(1,length(LoalNeu{1})),0.5*ones(1,length(LoalNeu{2})),0.8*ones(1,length(LoalNeu{3}))],seg_ind)
%     end
    switch i
        case 1
            ind = 11:21:num_files;
        case 2
            ind = 1:21:num_files;
        case 3
            ind = 21:21:num_files;
    end
    Fr = zeros(3,length(2.3e4:1e5));
    for j = 1:NumP
        MovVa = zeros(100,length(R.spike_hist{1}));
        for id_out = 1:100 % 100*(NumP-1)+1:100*NumP %  % 11:21:num_files
            %         fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
            fprintf('\t File name: %s\n', files{ind(id_out)});
            R = load(files{ind(id_out)});
            MovVa(id_out,:) = sum(movsum(full(R.spike_hist{1}(LoalNeu{j},:)),bin,2))/bin*1e4; %
        end
        Fr(j,:) = smooth(mean(MovVa(:,2.3e4:1e5),1),1e4)';
    end
%     subplot(3,3,3*i-1)
%     imagesc((2.3e4:length(Fr)-2e4)*1e-4,1,Fr(1,2.3e4:end-2e4))
%     hold on
%     imagesc((2.3e4:length(Fr)-2e4)*1e-4,3,Fr(2,2.3e4:end-2e4))
%     hold on
%     imagesc((2.3e4:length(Fr)-2e4)*1e-4,5,Fr(3,2.3e4:end-2e4))
%     ylim([0.5 5.5])
% %     c = gray;
% %     c = flipud(c);
% %     colormap(c);
%     colorbar
%     if i == 1
%         c1 = caxis;
%         title('Population Firing Rate(Hz)')
%     elseif i == 3
%         c2 = caxis;
%         cn = [min([c1 c2]), max([c1 c2])];
%         caxis(cn)
%         subplot(3,3,2)
%         caxis(cn)
%         subplot(3,3,5)
%         caxis(cn)
%         subplot(3,3,8)
%     end
%     xlabel('Time(s)','fontSize',10)
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
    subplot(3,3,3*i)
    plot((2.3e4:length(Fr)-2e4)*1e-4,Fr(2,2.3e4:end-2e4))
    xlabel('Time(s)','fontSize',10)
    ylabel('FR(Hz)','fontSize',10)
    if i == 2
        ylim([30 520])
    else
        ylim([30 100])
    end
end
%%
set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc WMFig1 % this is the trick!!
end