function LFPAmpPatternAnalysis2
dir_strut = dir('3DBurstLFPCut0*minTime30SR1000P95.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
Slope = [];
for i = 1:num_files
    R = load(files{i});
%     % Size of Bursts v.s. Time Interval between Bursts
%     PS = R.patternScale(1:end-1);
%     ESInterval = R.centInterval(2:end); % end-start interval:ms
% %     ESInterval = 5*ESInterval;
%     dat = [PS',ESInterval'];
%     n = hist3(dat,[50 50]);
%     n1 = n';
%     n1(size(n,1) + 1, size(n,2) + 1) = 0;
%     xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
%     yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
%     figure(1)
%     subplot(1,2,1)
%     h = pcolor(xb,yb,n1);
%     set(h, 'EdgeColor', 'none');
%     h.ZData = ones(size(n1)) * -max(max(n));
%     colormap(hot)
%     oldcmap = colormap(gray);
%     colormap( flipud(oldcmap) );
%     colorbar
%     [r,p] = corrcoef(PS',ESInterval');
%     if p(2) > 10^(-4)
%         pvalue = ['p = ',num2str(p(2),'%.4f')];
%     else
%         pvalue = ['p < 10^{-4}'];
%     end
%     tp = {['r = ',num2str(r(2),'%.4f')],pvalue};
%     text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
%     xlabel('Pattern Size')
%     ylabel('Interval among Bursts(ms)')
%     % Size of Bursts v.s. Bursts Jumping Distance
%     ESDist = 30*R.distCent(2:end); % end-start distance:um
%     dat = [PS',ESDist'];
%     n = hist3(dat,[50 50]);
%     n1 = n';
%     n1(size(n,1) + 1, size(n,2) + 1) = 0;
%     xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
%     yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
%     subplot(1,2,2)
%     h = pcolor(xb,yb,n1);
%     set(h, 'EdgeColor', 'none');
%     h.ZData = ones(size(n1)) * -max(max(n));
%     colormap(hot)
%     oldcmap = colormap(gray);
%     colormap( flipud(oldcmap) );
%     colorbar
%     [r,p] = corrcoef(PS',ESDist');
%     if p(2) > 10^(-4)
%         pvalue = ['p = ',num2str(p(2),'%.4f')];
%     else
%         pvalue = ['p < 10^{-4}'];
%     end
%     tp = {['r = ',num2str(r(2),'%.4f')],pvalue};
%     text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
%     xlabel('Pattern Size')
%     ylabel('Jumping Distance between Bursts(um)')
%     
%     figure(2)
%     % Distribution of Burst Size
%     subplot(2,2,1)
%     histogram(R.patternScale,20)
%     xlabel('Pattern Size')
%     ylabel('Count')
%     subplot(2,2,2)
%     [N,edges] = histcounts(R.patternScale,20);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Pattern Size')
%     ylabel('Count')
%     % Distribution of Burst Duration (over shreshold)
%     subplot(2,2,3)
%     histogram(R.Duration,20)
%     xlabel('Pattern Duration(ms)')
%     ylabel('Count')
%     subplot(2,2,4)
%     [N,edges] = histcounts(R.Duration,20);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Pattern Duration(ms)')
%     ylabel('Count')
%     
%     % Distribution of Propagation Step within Bursts
%     if ~isfield(R,'PStep')
%         PStep = [];
%         for j = 1:length(R.WCentroids)
%             centroids = R.WCentroids{j};
%             for k = 1:length(centroids)-1
%                 PStep = [PStep Distance_xy(centroids(k,2),centroids(k,3),centroids(k+1,2),centroids(k+1,3),20)];
%             end
%         end
%         save(files{i},'PStep','-append')
%     else
%         PStep = R.PStep;
%     end
%     figure(3)
%     subplot(3,2,1)
%     histogram(30*PStep,50) % um
%     xlabel('Moving Distance per Step(um)')
%     ylabel('Count')
%     subplot(3,2,2)
%     [N,edges] = histcounts(30*PStep,100);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Moving Distance per Step(um)')
%     ylabel('Count')
%     % Distribution of Burst Jumping Distance
%     subplot(3,2,3)
%     histogram(30*R.distCent(2:end),50)
%     xlabel('Jumping Distance(um)')
%     ylabel('Count')
%     subplot(3,2,4)
%     [N,edges] = histcounts(30*R.distCent(2:end),100);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Jumping Distance(um)')
%     ylabel('Count')
%     % Distribution of Propagation Step within Bursts + Jumping Distance between Bursts
%     AllMove = 30*[PStep R.distCent(2:end)]; % um
%     subplot(3,2,5)
%     histogram(AllMove,50)
%     xlabel('Moving Distance(um)')
%     ylabel('Count')
%     subplot(3,2,6)
%     [N,edges] = histcounts(AllMove,100);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Moving Distance(um)')
%     ylabel('Count')
%     
%     % Distribution of Burst Interval
%     figure(4)
%     subplot(1,2,1)
%     histogram(ESInterval,50)
%     xlabel('Burst Interval(ms)')
%     ylabel('Count')
%     subplot(1,2,2)
%     [N,edges] = histcounts(ESInterval,20);
%     Y = (edges(1:end-1)+edges(2:end))/2;
%     loglog(Y,N,'o')
%     xlabel('Burst Interval(ms)')
%     ylabel('Count')
    
    % MSD
    center = R.WCentroids ;
    numB = length(center);
%     figure(5)
    for iBurst = 1:numB % 20 % 1:numB
        Trajectory = [];
        TrajectoryTemp = center{iBurst}(:,[2:3 1]);
        Trajectory = [Trajectory;TrajectoryTemp] ;
        Trajectory(:,1:2) = (Trajectory(:,1:2) - 0.5)/63*600;
        [MSD,tau] = get_MSD_PBC(Trajectory(1:end,:));
%         tau = 5*tau;
        p = polyfit(log(tau(1:10)),log(MSD(1:10)),1);
        % select only 1 burst to plot
        if iBurst == 4 && i == 1
            subplot(1,2,1)
            loglog(tau,MSD,'.-')
            y = exp(polyval(p,log(tau))) ;
            hold on
            loglog(tau,y)
            title('Mean Square Distance versus \tau within a burst')
            xlabel('\tau (ms)')
            ylabel('Mean Square Distance (electrode)')
            str = {'p = ',num2str(p(1))};
            text(max(tau),max(y),str)
        end
        slope(iBurst) = p(1) ;
    end
    Slope = [Slope slope];
end

    subplot(1,2,2)
    histogram(Slope,20)
    disp(nanmean(Slope))
    xlabel('Slope')
    ylabel('Count')
    ts = sprintf('Mean Slope = %.4f ', nanmean(Slope));
    title(ts);
    %%
%     center = WCentroids ;
%     numB = length(center);
%     for iBurst = 1:numB % 20 % 1:numB
%         Trajectory = [];
%         TrajectoryTemp = center{iBurst}(:,[2:3 1]);
%         Trajectory = [Trajectory;TrajectoryTemp] ;
%         
%         [MSD,tau] = get_MSD(Trajectory(1:end,:));
%         p = polyfit(log(tau(1:10)),log(MSD(1:10)),1);
%         % select only 1 burst to plot
%         loglog(tau,MSD,'.-')
%         y = exp(polyval(p,log(tau))) ;
%         hold on
%         loglog(tau,y)
%         title('Mean Square Distance versus \tau within a burst')
%         xlabel('\tau (ms)')
%         ylabel('Mean Square Distance (electrode)')
%         str = {'p = ',num2str(p(1))};
%         text(max(tau),max(y),str)
%         next = input('\t Next figure?');
%         delete(gcf);
%     end
%     %%
%     % Pattern size autocorrelation
%     figure(6)
%     j = 15;
%     subplot(2,1,1)
%     plot(R.instantScale{j}(:,2))
%     xlabel('Time(ms)')
%     ylabel('Instantaneous Pattern Size')
%     subplot(2,1,2)
%     autocorr(R.instantScale{j}(:,2),length(R.instantScale{j})-1)
%     xlabel('Time(ms)')
%     ylabel('Autocorrelation')
% end
% %% summing all trial results in pattern size
% patternScale = [];
% for i = 1:num_files
%     R = load(files{i});
%     patternScale = [patternScale R.patternScale];
% end
% subplot(2,2,1)
% histogram(patternScale)
% xlabel('Pattern Size')
% ylabel('Count')
% subplot(2,2,2)
% [N,edges] = histcounts(patternScale,60);
% Y = (edges(1:end-1)+edges(2:end))/2;
% loglog(Y,N,'o')
% hold on;
% Y = Y(N > 0);
% N = N(N > 0);
% Y = Y(7:end);
% N = N(7:end);
% v = polyfit(log10(Y),log10(N),1);
% x = Y;
% y = 10^v(2)*x.^v(1);
% loglog(x,y)
% str = ['p = ',num2str(v(1))];
% text(min(x),max(y),str)
% xlabel('Pattern Size')
% ylabel('Count')
% % summing all trial results in pattern duration
% Duration = [];
% for i = 1:num_files
%     R = load(files{i});
%     Duration = [Duration R.Duration];
% end
% % Duration = 5*Duration;
% subplot(2,2,3)
% histogram(Duration)
% xlabel('Pattern Duration')
% ylabel('Count')
% subplot(2,2,4)
% [N,edges] = histcounts(Duration,60);
% Y = (edges(1:end-1)+edges(2:end))/2;
% loglog(Y,N,'o')
% hold on;
% Y = Y(N > 0);
% N = N(N > 0);
% % Y = Y(1:end-8);
% % N = N(1:end-8);
% v = polyfit(log10(Y),log10(N),1);
% x = Y;
% y = 10^v(2)*x.^v(1);
% loglog(x,y)
% str = ['p = ',num2str(v(1))];
% text(mean(x),max(y),str)
% xlabel('Pattern Duration')
% ylabel('Count')
% %%
% patternScale = [];
% for i = 1:num_files
%     R = load(files{i});
%     patternScale = [patternScale R.patternScale];
% end
% Duration = [];
% for i = 1:num_files
%     R = load(files{i});
%     Duration = [Duration R.Duration];
% end
% %%
% patternScale = R.patternScale;
% Duration = R.Duration;
% [alpha, xmin, L]=plfit(patternScale);
% plplot(patternScale, xmin, alpha)
% % [p,gof]=plpva(patternScale, xmin)
% str = ['SizeDistribution p = ',num2str(-alpha)];
% title(str)
% [alpha, xmin, L]=plfit(Duration);
% plplot(Duration, xmin, alpha)
% % [p,gof]=plpva(Duration, xmin)
% str = ['DurationDistribution p = ',num2str(-alpha)];
% title(str)
% %% Maximum likelihood estimation
% % power law distribution
% % a is the minimum value
% signal = patternScale; % (patternScale>2890);
% a = 2900 ;
% pf_power = @(x,alpha) x.^-alpha*(alpha-1)*a^(alpha-1);
% [lambdaHat,lambdaCI] = mle(signal, 'pdf',pf_power, 'start',[1.1], 'lowerbound',[1], 'upperbound', [3.5]) ;
% L = sum(log(pf_power(signal,lambdaHat(1)))) ;
% 
% % % exponential distribution
% pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*a) ;
% [lambdaHat3,lambdaCI] = mle(signal, 'pdf',pf_exp, 'start',[0] , 'lowerbound',[0], 'upperbound', [2]) ;
% L3 = sum(log(pf_exp(signal,lambdaHat3(1)))) ;
% 
% % signal = Duration; % (patternScale>3300);
% % a = 30 ;
% % pf_power = @(x,alpha) x.^-alpha*(alpha-1)*a^(alpha-1);
% % [lambdaHat,lambdaCI] = mle(signal, 'pdf',pf_power, 'start',[1.1], 'lowerbound',[1], 'upperbound', [7]) ;
% % L = sum(log(pf_power(signal,lambdaHat(1)))) ;
% % 
% % % % exponential distribution
% % pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*a) ;
% % [lambdaHat3,lambdaCI] = mle(signal, 'pdf',pf_exp, 'start',[0] , 'lowerbound',[0], 'upperbound', [2]) ;
% % L3 = sum(log(pf_exp(signal,lambdaHat3(1)))) ;
% 
% % [lambdaHat3,lambdaCI] = mle(signal,  'distribution','exponential') ;
% %% Maximum likelihood estimation 2 variables
% pf_power = @(x,alpha,a) x.^-alpha*(alpha-1)*a^(alpha-1);
% [lambdaHat,lambdaCI] = mle(patternScale, 'pdf',pf_power, 'start',[1.1,500], 'lowerbound',[1,10], 'upperbound', [5,1e4]) ;
% L = sum(log(pf_power(patternScale,lambdaHat(1), lambdaHat(2) ))) ;
% 
% % exponential distribution
% pf_exp = @(x,lambda,a) exp(-lambda*x) *lambda*exp(lambda*a) ;
% [lambdaHat3,lambdaCI] = mle(patternScale, 'pdf',pf_exp, 'start',[0,500] , 'lowerbound',[0,30], 'upperbound', [5,1e4]) ;
% L3 = sum(log(pf_exp(patternScale,lambdaHat3(1),lambdaHat3(2)))) ;
% 
% % gaussian distribution
% pf_gaussian = @(x,mu,sigma) exp(-lambda*x) *lambda*exp(lambda*a) ;
% [lambdaHat5,lambdaCI] = mle(patternScale, 'pdf',pf_exp, 'start',[0,500] , 'lowerbound',[0,30], 'upperbound', [5,1e4]) ;
% L3 = sum(log(pf_gaussian(patternScale,lambdaHat5(1),lambdaHat5(2)))) ;
% 
% % pf_power = @(x,alpha,a) x.^-alpha*(alpha-1)*a^(alpha-1);
% % [lambdaHat,lambdaCI] = mle(Duration, 'pdf',pf_power, 'start',[1.1,30], 'lowerbound',[1,10], 'upperbound', [5,100]) ;
% % L = sum(log(pf_power(Duration,lambdaHat(1), lambdaHat(2) ))) ;
% % 
% % % exponential distribution
% % pf_exp = @(x,lambda,a) exp(-lambda*x) *lambda*exp(lambda*a) ;
% % [lambdaHat3,lambdaCI] = mle(Duration, 'pdf',pf_exp, 'start',[0,30] , 'lowerbound',[0,30], 'upperbound', [5,100]) ;
% % L3 = sum(log(pf_exp(Duration,lambdaHat3(1),lambdaHat3(2)))) ;
% %% plot
% close all
% z = linspace(2000,9210,50); % patternScale:580,9210; Duration:30,100,50
% x = 0.5*(z(1:end-1)+z(2:end));
% subplot(2,2,1)
% histogram(signal,z,'Normalization','pdf')
% hold on
% plot(x,pf_power(x,lambdaHat))
% subplot(2,2,2)
% c = histcounts(signal,z,'Normalization','pdf');
% loglog(x,c,'.','markersize',20)
% hold on
% loglog(x,pf_power(x, lambdaHat)); 
% str = ['p = ',num2str(-lambdaHat),',L = ',num2str(L)];
% title(str)
% subplot(2,2,3)
% histogram(signal,z,'Normalization','pdf')
% hold on
% plot(x,pf_exp(x,lambdaHat3))
% xlabel('Size')
% subplot(2,2,4)
% c = histcounts(signal,z,'Normalization','pdf');
% loglog(x,c,'.','markersize',20)
% hold on
% loglog(x,pf_exp(x, lambdaHat3));
% xlabel('Size')
% str = ['L = ',num2str(L3)];
% title(str)
% %%
% % figure
% % nums = nums(duration > 10);
% % duration = duration(duration > 10);
% % subplot(2,2,1)
% % histogram(nums)
% % xlabel('Pattern Size')
% % ylabel('Count')
% % subplot(2,2,2)
% % [N,edges] = histcounts(nums,60); % 60
% % Y = (edges(1:end-1)+edges(2:end))/2;
% % loglog(Y,N,'o')
% % hold on;
% % Y = Y(N > 0);
% % N = N(N > 0);
% % Y = Y(3:end); % 2:end-20
% % N = N(3:end);
% % v = polyfit(log10(Y),log10(N),1);
% % x = Y;
% % y = 10^v(2)*x.^v(1);
% % loglog(x,y)
% % str = ['p = ',num2str(v(1))];
% % text(mean(x),mean(y),str)
% % xlabel('Pattern Size')
% % ylabel('Count')
% % subplot(2,2,3)
% % histogram(duration)
% % xlabel('Pattern Duration(ms)')
% % ylabel('Count')
% % subplot(2,2,4)
% % [N,edges] = histcounts(duration,60); % 60
% % Y = (edges(1:end-1)+edges(2:end))/2;
% % loglog(Y,N,'o')
% % hold on;
% % Y = Y(N > 0);
% % N = N(N > 0);
% % % Y = Y(1:end); % 2:end-1
% % % N = N(1:end);
% % v = polyfit(log10(Y),log10(N),1);
% % x = Y;
% % y = 10^v(2)*x.^v(1);
% % loglog(x,y)
% % str = ['p = ',num2str(v(1))];
% % text(mean(x),mean(y),str)
% % xlabel('Pattern Duration(ms)')
% % ylabel('Count')
% % %%
% % % [alpha, xmin, L]=plfit(nums);
% % % plplot(nums, xmin, alpha)
% % % str = ['SizeDistribution p = ',num2str(-alpha)];
% % % title(str)
% % [alpha, xmin, L]=plfit(duration,'xmin',10);
% % plplot(duration, xmin, alpha)
% % str = ['DurationDistribution p = ',num2str(-alpha)];
% % title(str)
% %% Visualizing Instantaneous LFP Amplitude
% load(files{2},'Centroids');
% timev = [];
% for i = 1:length(Centroids)
%     timev = [timev Centroids{i}(:,1)'];
% end
% load('0011-LocalSpikesGrid.mat','SpikesGrid');
% % gamma_power_grid = gamma_power_grid(:,:,1:10:end);
% gamma_power_grid = SpikesGrid;
% Y = prctile(gamma_power_grid(:),95);
% % vidObj = VideoWriter('th30fs1e3GammaBurstAmpOver95.avi');
% % vidObj.Quality = 100;
% % vidObj.FrameRate = 6; % number of frames to display per second
% % open(vidObj);
% for vt = (1:1e5)+1.01e3
%     A = gamma_power_grid(:,:,vt);
%     GGrid = zeros(10);
%     GGrid(A >= Y) = 1;
%     currentBurst = A.*GGrid;
%     if ~ismember(vt,timev) % This is for minimum length limit for bursts
%         currentBurst = zeros(10);
%     end
%     imagesc(currentBurst)
%     colorbar
%     ts = sprintf('time = %8.1f ms', vt);
%     title(ts);
%     pause(0.5);
% %     writeVideo(vidObj, getframe(gca));
% end
% % close(gcf);
% % close(vidObj);
% end