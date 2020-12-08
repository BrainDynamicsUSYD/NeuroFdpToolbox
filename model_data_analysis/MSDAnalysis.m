% function MSDAnalysis
% adapt from Xian's MSDnew.m function
% ATTENTION: manually modify the number of electrodes below
dir_strut = dir('3DBurst3000*minTime30SR1000.mat'); % 3DBurstLFP000*minTime30SR1000P95.mat % 3DBurst3000*minTime30SR1000.mat
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 0 ;
slopeAll = [] ;
% Trajectory = [];
for i = 3:num_files
    fprintf('Loading 3DBurst.mat file %s...\n', files{i});
    R = load(files{i});
    center = R.WCentroids;
%     a = find(R.Duration>100);
    for iBurst = 1:size(center,2)
        Trajectory = [center{iBurst}(:,[2:3 1])] ;
        Trajectory(:,1:2) = Trajectory(:,1:2)/40*600; %% manually modify here according to electrodes!!! %%
        [MSD,tau] = get_MSD_PBC(Trajectory);
        tempMSD = [MSD,tau] ;
        p_temp = [];
        norErr = [] ;
        maxStart = 20 ;               % minimum tau used to fit MSD
        tau_max = maxStart:min(40,size(MSD,1)) ;
        for fitIdx = 1:length(tau_max)
            fitRange = (1:tau_max(fitIdx)) ;
            [pAll,~] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
            y = exp(polyval(pAll,log(tau(fitRange)))) ;
            errorRate = mean(abs((MSD(fitRange)-y(fitRange))./MSD(fitRange))) ;
            p_temp(fitIdx,:) = pAll ;
            norErr(fitIdx) = errorRate ;
        end
        %         % discard the MSD with error rate more than 10%
        if (min(norErr)<0.1) 
            bestIdx = find(norErr == min(norErr)) ;
            pBest = p_temp(bestIdx(end),:) ;
            %             y = exp(polyval(pBest,log(tau))) ;
            slopeAll = [slopeAll pBest(1)];
            if pBest(1) > 1.33 && pBest(1) < 1.34  % iBurst == 9 && i == 1
            figure
%             subplot(1,2,1)
            loglog(tau,MSD,'o')
            hold on
            loglog(tau(fitRange),y,'LineWidth',1.5)% *2)
            title('Mean Square Distance versus \tau') % within a burst')
            xlabel('\tau (ms)')
            ylabel('Mean Square Distance (um^2)')
            str = ['p = ',num2str(pBest(1))];
            text(max(tau),max(y),str)
%             subplot(1,2,2)
%             Trajectory = [center{iBurst}(:,2:3)] ;
%             plot(Trajectory(:,1),Trajectory(:,2),'.-')
%             xlim([0 40])
%             ylim([0 40])
%             title('Trajectory')
            next = input('\t Next figure?');
            delete(gcf);
            end
        else
            continue
        end
        %         % keep the superdiffusive one to do average
        if pBest(1)>1.1
            msdSuper = [msdSuper;tempMSD] ;
            countSuper = countSuper + 1;
        end
    end
end
%% Spikes pattern version
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
meanMSD = [] ;
stdMSD = [] ;
msdSuper = [] ;
countSuper = 0 ;
slopeAll = [] ;
hw = 31;
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    R = get_grid_firing_centre(R);
    R = {R};
    SaveRYG(R);
    disp('Done');
    R = R{1};
    Ind = find(~isnan(R.grid.quick.radius));
    t_mid = R.grid.t_mid(Ind);
    centre = R.grid.quick.centre(:,Ind);
    div = 100; % ~ms
    Spikespattern = cell(1,ceil(length(Ind)/div));
    for i = 1:length(Ind)
        Spikespattern{ceil(i/div)} = [Spikespattern{ceil(i/div)};t_mid(i) centre(:,i)'];
    end
    for iPattern = 1:length(Spikespattern)
        Trajectory = [Spikespattern{iPattern}(1:1:end,[2:3 1])] ;
        Trajectory(:,1:2) = Trajectory(:,1:2)/63*600; %% manually modify here according to electrodes!!! %%
        Trajectory(:,3) = 0.1*Trajectory(:,3);
        [MSD,tau] = get_MSD_PBC(Trajectory);
        tempMSD = [MSD,tau] ;
        p_temp = [];
        norErr = [] ;
        maxStart = 20 ;               % minimum tau used to fit MSD
        tau_max = maxStart:min(40,size(MSD,1)) ;
        for fitIdx = 1:length(tau_max)
            fitRange = (1:tau_max(fitIdx)) ;
            [pAll,~] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
            y = exp(polyval(pAll,log(tau(fitRange)))) ;
            errorRate = mean(abs((MSD(fitRange)-y(fitRange))./MSD(fitRange))) ;
            p_temp(fitIdx,:) = pAll ;
            norErr(fitIdx) = errorRate ;
        end
        if (min(norErr)<0.1)
            bestIdx = find(norErr == min(norErr)) ;
            pBest = p_temp(bestIdx(end),:) ;
            y = exp(polyval(pBest,log(tau))) ;
            slopeAll = [slopeAll pBest(1)];
%             figure
%             subplot(1,2,1)
%             loglog(tau,MSD,'.-')
%             hold on
%             loglog(tau(fitRange),y(fitRange))% *2)
%             title('Mean Square Distance versus \tau') % within a burst')
%             xlabel('\tau (ms)')
%             ylabel('Mean Square Distance (um^2)')
%             str = {'p = ',num2str(pBest(1))};
%             text(max(tau(fitRange)),max(y(fitRange)),str)
%             subplot(1,2,2)
%             Trajectory = [Spikespattern{iPattern}(:,2:3)] ;
%             plot(Trajectory(:,1),Trajectory(:,2),'.-')
%             xlim([-hw hw])
%             ylim([-hw hw])
%             title('Trajectory')
%             %         disp(iPattern)
%             next = input('\t Next figure?');
%             delete(gcf);
        else
            continue
        end
        %         % keep the superdiffusive one to do average
        if pBest(1)>1.1
            msdSuper = [msdSuper;tempMSD] ;
            countSuper = countSuper + 1;
        end
    end
end
%% plot for average MSD and distribution
% figure_width = 8.4; %cm
% figure_hight = 8.4; %cm
% figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', ...
%     'color','w', 'position', [0, 0, figure_width, figure_hight], ...
%     'PaperSize', [figure_width, figure_hight]); % this is the trick!
figure;
% subplot(1,2,1)
sigIn = msdSuper ;
[MSDsort,~,MSDIdx] = unique(sigIn(:,2)) ;
for idx = 1:length(MSDsort)
    meanMSD(idx) = mean(sigIn(find(MSDIdx==idx),1)) ;
    %     stdMSD(idx) = std(sigIn(find(MSDIdx==idx),1)) ;
end
plot(MSDsort,meanMSD,'o')
% errorbar(MSDsort,meanMSD,stdMSD)
set(gca,'yscale','log','xscale','log')
hold on
p = polyfit(log(MSDsort(1:20)'),log(meanMSD(1:20)),1) ;
tau = MSDsort(1:50);
y = exp(polyval(p,log(tau))) ;
loglog(tau,y,'LineWidth',1.5)
title(['Mean Diffusive MSD for SuperD (',num2str(countSuper),')'])
% set(gca,'FontSize', 15);
xlabel('\tau (ms)')
ylabel('Mean Square Distance (um^2)')
str = ['p = ',num2str(p(1))];
text(max(tau),max(y),str)
% text(-0.2,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% subplot(1,2,2)
% hist(slopeAll,40)
% xlabel('MSD Slope')
% ylabel('Count')
% ts = sprintf('Mean Slope = %.4f ', nanmean(slopeAll));
% title(ts);
% text(-0.2,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')

% set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
% print -depsc figure_size_control % this is the trick!! 
% end
