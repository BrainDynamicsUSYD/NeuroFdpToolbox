% Gamma Pattern Paper Figure 3
figure_width = 17.8; % cm
figure_hight = 9; % cm
figure('NumberTitle','off','name', 'GammaPaperFig3V3', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!
%%
load('0015-201806291103-38833_in_1530234329865_0_neurosamp.mat','LFP_grid');
R = load('0015-201806291103-38833_in_1530234329865_out_RYG.mat');
load('3DBurst30015minTime30SR1000.mat','WCentroids')
RR.LFP.LFP{1} = reshape(LFP_grid,63^2,[]);
[RR] = GetBurst2(RR);
sigBinary = reshape(RR.LFP.GammaBurstEvent.is_burst,63,63,[]);
clear LFP_grid
LFP_grid = reshape(RR.LFP.LFP_gamma_hilbert_abs,63,63,[]);
clear RR
S = full(R.spike_hist{1});
S = reshape(S,63,63,[]);
%%
pick = [0;2];% [35 0;20 40]; % [25 65;35 10]/[25 0;35 80];
j = 1;
for i = 1 % :13 % [1 9] % [1 4]
%     subplot(2,4,j)
    t = WCentroids{i}(end-pick(2,j),1)/R.dt;
    LFP = LFP_grid(:,:,t);
    clim = minmax(reshape(LFP_grid(:,:,(WCentroids{i}(1,1):WCentroids{i}(end,1))/R.dt),1,[]));
    imagesc(LFP'); % ,clim);
%     colorbar
    hold on
    [Row,Col] = find(sigBinary(:,:,t));
    plot(Row,Col,'k.','MarkerSize',1) %,'color',[1 0.8 0.8])
    hold on
    for tms = WCentroids{i}(1,1)+pick(1,j):WCentroids{i}(end,1)-pick(2,j)
        if tms == WCentroids{i}(1,1)+pick(1,j)
            temp = [63-WCentroids{i}(1+pick(1,j),2),WCentroids{i}(1+pick(1,j),3)]; % 63-WCentroids{i}(1,2:3);
            ind = 1+pick(1,j);
        else
            ind = ind + 1;
            y = [temp(2) WCentroids{i}(ind,3)];
            x = [temp(1) 63-WCentroids{i}(ind,2)];
            hold on
            plot(x,y,'r-','LineWidth',1) % ,'MarkerSize',5)
            temp = [63-WCentroids{i}(ind,2),WCentroids{i}(ind,3)]; % 63-WCentroids{i}(ind,2:3);
            if tms == WCentroids{i}(end,1)-pick(2,j)
                plot(temp(1),temp(2),'r.','MarkerSize',10)
            end
        end
        xlim([0 63])
        ylim([0 63])
    end
    if j == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',12)
    end
%     j = j + 1;
% next = input('\t Next figure?');
% close all
end
%%
pick = [25 0;35 5]; % [25 65;35 10]/[25 0;35 80];
j = 1;
for i = [1 7] % [1 4]
    subplot(2,4,j)
    t = WCentroids{i}(end-pick(2,j),1)/R.dt;
    LFP = LFP_grid(:,:,t);
%     Spikes = sum(S(:,:,(t-25:t+24)+2e4),3);
%     [row,col] = find(Spikes);
    clim = minmax(reshape(LFP_grid(:,:,(WCentroids{i}(1,1):WCentroids{i}(end,1))/R.dt),1,[]));
    imagesc(LFP'); % ,clim);
    colorbar
    hold on
    [Row,Col] = find(sigBinary(:,:,t));
    plot(Row,Col,'k.','MarkerSize',1) %,'color',[1 0.8 0.8])
    hold on
%     plot(row,col,'w.','MarkerSize',1)
%     hold on;
    for tms = WCentroids{i}(1,1)+pick(1,j):WCentroids{i}(end,1)-pick(2,j)
        if tms == WCentroids{i}(1,1)+pick(1,j)
            temp = [63-WCentroids{i}(1+pick(1,j),2),WCentroids{i}(1+pick(1,j),3)]; % 63-WCentroids{i}(1,2:3);
            ind = 1+pick(1,j);
        else
            ind = ind + 1;
            y = [temp(2) WCentroids{i}(ind,3)];
            x = [temp(1) 63-WCentroids{i}(ind,2)];
            hold on
            plot(x,y,'r--','LineWidth',1) % ,'MarkerSize',5)
            temp = [63-WCentroids{i}(ind,2),WCentroids{i}(ind,3)]; % 63-WCentroids{i}(ind,2:3);
            if tms == WCentroids{i}(end,1)-pick(2,j)
                plot(temp(1),temp(2),'r.','MarkerSize',10)
            end
        end
        xlim([0 63])
        ylim([0 63])
    end
    if j == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',12)
    end
    j = j + 1;
end
%%
cd ../GammaOscillation1600Electrodes/
dir_strut = dir('3DBurst3000*minTime30SR1000.mat'); 
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
for i = 1:num_files
    fprintf('Loading 3DBurst.mat file %s...\n', files{i});
    R = load(files{i});
    center = R.WCentroids;
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
        % discard the MSD with error rate more than 10%
        if (min(norErr)<0.1)
            bestIdx = find(norErr == min(norErr)) ;
            pBest = p_temp(bestIdx(end),:) ;
            slopeAll = [slopeAll pBest(1)];
            if i == 3 && iBurst == 296  % iBurst == 9 && i == 1
                subplot(2,3,4)
                loglog(tau,MSD,'.','MarkerSize',10)
                hold on
                loglog([tau(fitRange);100],[y;exp(polyval(pAll,log(100)))],'LineWidth',2)% *2)
                xlabel('\tau (ms)','fontsize',8)
                ylabel('MSD (um^2)','fontsize',8)
                str = ['p = ',num2str(pBest(1))];
                text(5,max(y),str,'fontsize',8)
                text(-0.2,1.02,'C','Units', 'Normalized','FontSize',12)
            end
        else
            continue
        end
        % keep the superdiffusive one to do average
        if pBest(1)>1.1
            msdSuper = [msdSuper;tempMSD] ;
            countSuper = countSuper + 1;
        end
    end
end

subplot(2,3,5)
hist(slopeAll,15)
xlabel('MSD exponent','fontsize',8)
ylabel('Count','fontsize',8)
ts = sprintf('Mean Alpha = %.2f ', nanmean(slopeAll));
% title(ts);
text(-0.2,1.02,'D','Units', 'Normalized','FontSize',12)
%%
hw = 31;
Dis = [];
eta = 0.6; % 0.9 % 0.67
deltaTime = [2 4 8 16 32]; % 30 ;
Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
ind = 1;
% hold on
subplot(4,2,[7,8])
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            DisplaceTemp = sqrt(sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2));
            Displace = [Displace;DisplaceTemp] ;
        end
    end
%     Displace = sqrt(Displace);
    nEdge = logspace(min(Displace),max(Displace),20) ;
    [n,nEdge] = histcounts(Displace,nEdge,'normalization','probability');
    loglog((nEdge(1:end-1)+nEdge(2:end))/2,n/sum(n),'.','MarkerSize',6,'color',Color(ind,:))
    hold on;
%     pd = fitdist(Displace,'Stable')
%     y = pdf(pd,nEdge);
%     loglog(nEdge,y/sum(y),'-','LineWidth',1,'color',Color(ind,:));
%     hold on ;
%     Dis = [Dis; Displace];
    ind = ind + 1;
end
xlabel('Displacement d(\tau) (\mum)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d(\tau))','fontsize',8)
legend({'\tau = 2 ms','\tau = 4 ms','\tau = 8 ms','\tau = 16 ms','\tau = 32 ms',},'fontsize',8,'edgecolor','none','color','none')
text(-0.2,1.02,'E','Units', 'Normalized','FontSize',12)
%%
ax=axes('Position',[0.85 0.35 0.1 0.1],'Unit','normalize',...
    'parent',1);
box on
Dis = [];
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            DisplaceTemp = sqrt(sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2));
            Displace = [Displace;DisplaceTemp] ;
        end
    end
    rmsd = deltaT^eta;
%     Displace = sqrt(Displace);
    Displace = Displace/rmsd;
    nEdge = logspace(log10(min(Displace)),log10(max(Displace)),20) ;
    [n,nEdge] = histcounts(Displace,nEdge,'normalization','probability');
    loglog((nEdge(1:end-1)+nEdge(2:end))/2,n/sum(n),'.','MarkerSize',6)
    hold on;
    Dis = [Dis; Displace];
end
% x2 = [linspace(-200,-80,100) x2 linspace(80,200,100)];
Dis = [-Dis;Dis];
% Dis3 = [-Dis3;Dis3];
% Dis2 = [-Dis2;Dis2];
x1 = linspace(0,0.35*max(Dis),200);
x2 = linspace(0,max(Dis),200);
pd = fitdist(Dis,'normal') % normal
y = pdf(pd,x1); % x1
pd = fitdist(Dis,'Stable')
yy = pdf(pd,x2);
loglog(x1,y/sum(y),'k--','LineWidth',1);
hold on ;
loglog(x2,yy/sum(yy),'r','LineWidth',2);
legend({'\tau = 2ms','\tau = 4 ms','\tau = 8 ms','\tau = 16 ms','\tau = 32 ms','Gaussian','\alpha stable'},'fontsize',8,'edgecolor','none','color','none')
xlabel('d_s(\tau)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d_s(\tau))','fontsize',8)
% text(-0.2,1.02,'E','Units', 'Normalized','FontSize',12)
%%
dir_strut2 = dir('3DBurst3000*minTime0SR1000.mat'); % 3DBurstLFPCut0*minTime30SR1000P95.mat; 3DBurst20*minTime30SR1000.mat
num_files = length(dir_strut2);
files2 = cell(1,num_files);
for id_out = 1:num_files
    files2{id_out} = dir_strut2(id_out).name;
end
Duration = [];
Interval = [];
for i = 1:num_files
    fprintf('Loading out_RYG.mat file %s...\n', files2{i});
    R = load(files2{i});
    Duration = [Duration R.Duration]; % R.centInterval(2:end)]; % R.distCent(2:end)]; % R.Duration]; % R.patternScale;
    Interval = [Interval R.centInterval(2:end)]; % R.centInterval(2:end)]; % R.distCent(2:end)]; % R.Duration]; % R.patternScale;
end
%%
subplot(2,4,3)
Signal = Duration; % Interval(Interval > 0); % [-Interval,Interval]; % Duration; % [-Duration,Duration]
% numPts = 20;
% nEdge = logspace(log10(min(Signal)),log10(max(Signal)),numPts) ;
% [N,edges] = histcounts(Signal,nEdge,'normalization','pdf') ;
% Y = 10.^((log10(edges(1:end-1))+log10(edges(2:end)))/2);
% Signal = [-Signal,Signal];
% Y = Y(N > 0);
% N = N(N > 0);
% YY = Y(12:end);
% NN = N(12:end);
% v = polyfit(log10(YY),log10(NN),1);
% x = Y;
% y = 10^v(2)*x.^v(1);
% plot(Y,N,'.','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
% hold on
% plot(x,y,'LineWidth',2)
% legend({'Data','Power-law'},'fontsize',8)
[N,edges] = histcounts(Signal,50,'normalization','probability');
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
Signal = [-Signal,Signal];
pd = fitdist(Signal','Stable')
y = pdf(pd,Y);% 2*  *(Y(2)-Y(1))
hold on
loglog(Y,y/sum(y),'r','LineWidth',2);
pd2 = fitdist(Signal','normal')
y2 = pdf(pd2,Y);% 2*
hold on
loglog(Y,y2/sum(y2),'k--','LineWidth',2);
legend({'Data','Levy stable','Gaussian'},'fontsize',8)
xlim([0 4e2])
ylim([6e-5 0.2])
xlabel('Duration','fontsize',8)
ylabel('Probability {\itP} ({\itT})','fontsize',8)
text(-0.2,1.02,'B','Units', 'Normalized','FontSize',12)
% ax=axes('Position',[0.55 0.8 0.1 0.1],'Unit','normalize',...
%     'parent',1);
% box on
% loglog(Y,N,'.','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
% hold on
% loglog(Y,y/sum(y),'r','LineWidth',2);
% xlim([10 4e2])
% ylim([6e-6 1])
% set(gca,'XTickLabel',[],'YTickLabel',[])

%%
% subplot(1,2,1)
Signal = Duration(Duration > 0); % Interval(Interval > 0); % [-Interval,Interval]; % Duration; % [-Duration,Duration]
numPts = 20;
nEdge = logspace(log10(min(Signal)),log10(max(Signal)),numPts) ;
[N,edges] = histcounts(Signal,nEdge,'normalization','pdf') ;
Y = 10.^((log10(edges(1:end-1))+log10(edges(2:end)))/2);
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'.','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
% Signal = [-Signal,Signal];
pd = fitdist(Signal','Stable')
y = pdf(pd,Y);% 2*
hold on
loglog(Y,y,'r','LineWidth',2);
pd = fitdist(Signal','normal')
y = pdf(pd,Y);% 2*
hold on
loglog(Y,y,'k--','LineWidth',2);
legend({'Data','Stable','Gaussian'},'fontsize',8)
xlim([0 4e2])
ylim([6e-6 1])
xlabel('Duration T (ms)','fontsize',8)
ylabel('PDF','fontsize',8)
% text(-0.2,1.02,'B','Units', 'Normalized','FontSize',12)
%%
subplot(2,4,4)
Signal = Interval(Interval > 0); % Interval(Interval > 0); % [-Interval,Interval]; % Duration; % [-Duration,Duration]
% numPts = 20;
% nEdge = logspace(log10(min(Signal)),log10(max(Signal)),numPts) ;
% [N,edges] = histcounts(Signal,nEdge,'normalization','pdf') ;
% Y = 10.^((log10(edges(1:end-1))+log10(edges(2:end)))/2);
% % Signal = [-Signal,Signal];
% Y = Y(N > 0);
% N = N(N > 0);
% YY = Y(12:end);
% NN = N(12:end);
% v = polyfit(log10(YY),log10(NN),1);
% x = Y;
% y = 10^v(2)*x.^v(1);
% plot(Y,N,'.','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
% hold on
% plot(x,y,'LineWidth',2)
% legend({'Data','Power-law'},'fontsize',8)
[N,edges] = histcounts(Signal,100,'normalization','probability');
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
Signal = [-Signal,Signal];
pd = fitdist(Signal','Stable')
y = pdf(pd,Y);% 2* *(Y(2)-Y(1))
hold on
loglog(Y,y/sum(y),'r','LineWidth',2);
pd2 = fitdist(Signal','normal')
y2 = pdf(pd2,Y);% 2*
hold on
loglog(Y,y2/sum(y2),'k--','LineWidth',2);
legend({'Data','Levy stable','Gaussian'},'fontsize',8,'edgecolor','none','color','none')
xlim([0 7e2])
ylim([6e-5 0.1])
xlabel('Interval','fontsize',8)
ylabel('Probability {\itP} (\delta)','fontsize',8)
% text(-0.2,1.02,'B','Units', 'Normalized','FontSize',12)
% ax=axes('Position',[0.85 0.8 0.1 0.1],'Unit','normalize',...
%     'parent',1);
% box on
% loglog(Y,N,'.','MarkerSize',12) % semilogy(Y,N/sum(N),'.','MarkerSize',10)
% hold on
% loglog(Y,y/sum(y),'r','LineWidth',2);
% xlim([10 4e2])
% ylim([6e-6 1])
% set(gca,'XTickLabel',[],'YTickLabel',[])
%%
subplot(1,2,2)
Signal = Interval(Interval > 0); % Interval(Interval > 0); % [-Interval,Interval]; % Duration; % [-Duration,Duration]
nEdge = logspace(log10(min(Signal)),log10(max(Signal)),numPts) ;
[N,edges] = histcounts(Signal,nEdge,'normalization','pdf') ;
Y = 10.^((log10(edges(1:end-1))+log10(edges(2:end)))/2);
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'.','MarkerSize',12)
% Signal = [-Signal,Signal];
pd = fitdist(Signal','Stable')
y = pdf(pd,Y);% 2*
hold on
loglog(Y,y,'r','LineWidth',2);
pd = fitdist(Signal','normal')
y = pdf(pd,Y);% 2*
hold on
loglog(Y,y,'k--','LineWidth',2);
legend({'Data','Stable','Gaussian'},'fontsize',8)
xlim([0 7e2])
ylim([6e-6 1])
xlabel('Interval (ms)','fontsize',8)
ylabel('PDF','fontsize',8)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc GammaPaperFig3V3 % this is the trick!!