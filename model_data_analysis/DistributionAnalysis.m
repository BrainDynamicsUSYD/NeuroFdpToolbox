% function DistributionAnalysis
% analyze distribution: power law OR lognormal
% dir_strut = dir('*_out_RYG.mat');
dir_strut2 = dir('3DBurst3000*minTime0SR1000.mat'); % 3DBurstLFPCut0*minTime30SR1000P95.mat; 3DBurst20*minTime30SR1000.mat
num_files = length(dir_strut2);
% files = cell(1,num_files);
files2 = cell(1,num_files);
for id_out = 1:num_files
%     files{id_out} = dir_strut(id_out).name;
    files2{id_out} = dir_strut2(id_out).name;
end
% stx = 'Amplitude(a.u.)';
% sty = 'Count';
% stt = 'Amplitude Distribution';
%%
Duration = [];
Interval = [];
Size = [];
for i = 1:num_files
        fprintf('Loading out_RYG.mat file %s...\n', files2{i});
        R = load(files2{i});
%         R = GetGamma(R);
%         R = GetBurst(R);
    %     load(files2{i});
    %%% GET gamma cycle Amplitude
    %     R = GetICI(R);
    %     LFP_gamma = R.LFP.LFP_gamma;
    %     [no,~] = size(LFP_gamma); % R.LFP.LFP_broad; R.pop_stats.V_mean{1}
    %     % set(gcf,'color','w');
    %     for i = 1:no
    %         [p1,l1] = findpeaks(LFP_gamma(i,:));
    %         [p2,l2] = findpeaks(-LFP_gamma(i,:));
    %         l = min(length(l1),length(l2));
    %         m = 1;
    %         for j = 1:l
    %             while l1(j) >= l2(m) && m <= l
    %                 m = m + 1;
    %             end
    %             if m > l
    %                 break
    %             end
    %             Amp = [Amp p1(j) + p2(m)];
    %             m = m + 1;
    %         end
    %     end
    %%% Get single electrode burst Duration/Interval
%     for j = 1:length(R.LFP.GammaBurstEvent.flat_du_steps)
%         Duration = [Duration 0.1*R.LFP.GammaBurstEvent.burst_du_steps{j}]; % ms
%         Interval = [Interval 0.1*R.LFP.GammaBurstEvent.flat_du_steps{j}]; % ms
%         % Bursts & Spikes
% %         for k = 1:length(R.LFP.GammaBurstEvent.burst_start_steps{j})
% %             start = R.LFP.GammaBurstEvent.burst_start_steps{j}(k);
% %             d = R.LFP.GammaBurstEvent.burst_du_steps{j}(k);
% %             Spikes = sum(double(full(R.spike_hist{1}(:,start:start+d-1))));
% %             Duration = [Duration Spikes];
% %         end
%     end
    %%% Get Burst drift Frequency
%     for j = 1:length(Start)
%         Interval = [Interval Start{j}(2:end)-End{j}(1:end-1)]; % s
%     end
%     Interval = 1e3*Interval; % ms
%     Duration = Duration*1e3;
%     Duration(Duration < 15) = [];
%     Drift = [Dfre{:}];
    %%% Get burst only Amplitude
%     s = size(R.LFP.LFP{1});
%     AmpHil = zeros(s);
%     for j = 1:s(1)
%         AmpHil(j,:) = abs(hilbert(R.LFP.LFP_gamma(j,:)));
%     end
%     Y = prctile(AmpHil(:),95);
%     AmpHil(AmpHil < Y) = 0;
%     AmpHil(AmpHil > 0) = 1;
%     for j = 1:s(1)
%          temp = AmpHil(j,:);
%          [~,duration,interval,~,~] = seq_postprocess(temp,1,1);
%          Duration = [Duration 0.1*duration];
%          Interval = [Interval 0.1*interval];
%     end    
%     Duration(Duration < 15) = [];
    %%% Get spatiotemporal pattern Duration/Interval/Size
        Duration = [Duration R.Duration]; % R.centInterval(2:end)]; % R.distCent(2:end)]; % R.Duration]; % R.patternScale;
%         Duration = Duration(Duration > 0);
        Size = [Size R.patternScale]; % R.centInterval(2:end)]; % R.distCent(2:end)]; % R.Duration]; % R.patternScale;
%         Size = Size(Size > 0);
        Interval = [Interval R.centInterval(2:end)]; % R.centInterval(2:end)]; % R.distCent(2:end)]; % R.Duration]; % R.patternScale;
%         Interval = Interval(Interval > 0);
    %%% Get spatiotemporal burst Width-like:instantScale
    %     for j = 1:length(R.Width) % instantScale/Width
    %         width = [width R.Width{j}]; % R.instantScale{j}(:,2)'];
    %     end
    %%% Get sum of individual displacement in each burst
    %     for j = 1:length(R.WCentroids) % instantScale/Width
    %         %         SumD = 0;
    %         %         for k = 1:length(R.WCentroids{j})-1
    %         % %             SumD = SumD + Distance_xy(R.WCentroids{j}(k,2),R.WCentroids{j}(k,3),R.WCentroids{j}(k+1,2),R.WCentroids{j}(k+1,3),40);
    %         %             Duration = [Duration Distance_xy(R.WCentroids{j}(k,2),R.WCentroids{j}(k,3),R.WCentroids{j}(k+1,2),R.WCentroids{j}(k+1,3),40)];
    %         %         end
    %         %         SumD = SumD/(length(R.WCentroids{j})-1);
    %         %         Duration = [Duration R.distCent(2:end)]; %./R.centInterval(2:end)]; % SumD]; % R.instantScale{j}(:,2)'];
    %         %         Duration(Duration==Inf) = [];
    %         if j == 1
    %             pre = AveragePBC(R.WCentroids{j}(:,2:3),40);
    %             continue
    %         end
    %         next = AveragePBC(R.WCentroids{j}(:,2:3),40);
    %         Duration = [Duration Distance_xy(pre(1),pre(2),next(1),next(2),40)];
    %         pre = next;
    %     end
end
%%
    subplot(1,2,1)
    histogram(Duration)
    xlabel('Duration(ms)')
    ylabel('Count')
    subplot(1,2,2)
    histogram(Interval)
    xlabel('Interval(ms)')
    ylabel('Count')
%     subplot(1,3,3)
%     histogram(Drift)
%     xlabel('Drift(Hz)')
%     ylabel('Count')
%%
Y = prctile(gamma_power_grid(:),95);
s = size(gamma_power_grid);
mark = 0;
D = []; % propagate step size
L = []; % burst last time
xwidth = [];
ywidth = [];
finish = 1;
% subplot(1,3,1)
for t = 1:s(3)
    A = gamma_power_grid(:,:,t);
    [peak_mag,I] = max(A(:));
    [I_row, I_col] = ind2sub(s([1 2]),I);
    A_c = circshift(A, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
    GGrid = zeros(s([1 2]));
    GGrid(A_c >= Y) = 1;
    CC = bwconncomp(GGrid);
    validRegionIdx = 0;
    count = 1;
    for iRegion = 1:size(CC.PixelIdxList,2)
        if(size(CC.PixelIdxList{iRegion},1)>40)
            validRegionIdx(count) = iRegion ;
            count = count + 1 ;
        end
    end
    if count == 2
        B = regionprops(CC,'BoundingBox'); % Smallest rectangle containing the region
        Boundary = cat(1, B.BoundingBox);
        xwidth = [xwidth Boundary(3)];
        ywidth = [ywidth Boundary(4)];
    end
end
xwidth = 600/63*xwidth;  % um
ywidth = 600/63*ywidth;  % um

Duration(Duration <= 0) = [];
%%
Signal = Duration; % /40*600; % (xwidth+ywidth)/2;
subplot(1,2,1)
[N,edges] = histcounts(Signal,100);
Y = (edges(1:end-1)+edges(2:end))/2;
semilogx(Y,N/sum(N),'.')
parmhat = lognfit(Signal);
xlabel('Duration(ms)')
ylabel('Probability')
title('(semilog)')
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(Y,YY/sum(YY),'r')
legend('data','fit lognormal')
% text(-0.28,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% subplot(1,3,2)
% Signal = Size;
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% semilogx(Y,N/sum(N),'.')
% parmhat = lognfit(Signal);
% xlabel('Size(a.u.)')
% ylabel('Probability')
% title('(semilog)')
% hold on
% YY = lognpdf(Y,parmhat(1),parmhat(2));
% semilogx(Y,YY/sum(YY),'r')
% legend('data','fit lognormal')
subplot(1,2,2)
Signal = Interval; % /40*600; % (xwidth+ywidth)/2;
[N,edges] = histcounts(Signal,100);
Y = (edges(1:end-1)+edges(2:end))/2;
semilogx(Y,N/sum(N),'.')
parmhat = lognfit(Signal);
xlabel('Interval(ms)')
ylabel('Probability')
title('(semilog)')
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(Y,YY/sum(YY),'r')
legend('data','fit lognormal')
% text(-0.28,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%% NCC TOOLBOX
Signal = I_exc; % width(1:round(end/2));
[tau, xmin, xmax, L] = plmle(Signal,'xmin',1); % ,'xmax',16);
fitParams = struct;
% tau = 8;
fitParams.tau = tau;
fitParams.xmin = xmin;
fitParams.xmax = xmax;
plotParams = struct;
plotParams.dot = 'on';
titleString = ['pl:alpha = ',num2str(tau,'%.4f')];
% data = plplottool(Signal,'plotParams',plotParams);
data = plplottool(Signal,'plotParams',plotParams,'fitParams',fitParams, 'title',titleString);
%% Normal power law fit
% subplot(1,2,1)
% Signal = Duration;
figure
[N,edges] = histcounts(Signal,50);
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o')
hold on;
Y = Y(N > 0);
N = N(N > 0);
YY = Y(5:40);
NN = N(5:40);
v = polyfit(log10(YY),log10(NN),1);
x = Y(2:end);
y = 10^v(2)*x.^v(1);
loglog(x,y,'LineWidth',1.5)
str = ['\alpha = ',num2str(-v(1))];
text(min(x),2*max(y),str)
% hold on;
% Ym = Y(13:end);
% Nm = N(13:end);
% v = polyfit(log10(Ym),log10(Nm),1);
% x = Y(13:end);
% y = 10^v(2)*x.^v(1);
% loglog(x,y,'LineWidth',1.5)
% str = ['p = ',num2str(v(1))];
% text(min(x),2*max(y),str)
xlabel('Size') % Interval(ms) % Duration(ms) % Size
ylabel('Count')
%%
subplot(1,3,2)
Signal = Size;
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o')
xlabel('Size(a.u.)')
ylabel('Count')
subplot(1,3,3)
Signal = Interval;
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o')
xlabel('Interval(ms)')
ylabel('Count')
%%
powerIndx = 1:0.01:5 ;
x = 10.^powerIndx ;
[n, xout] = hist(Signal,x);
%[n, xout] = hist(patternScale,1000);
n=n/sum(n) ;
plot(xout,n,'.','MarkerSize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
title('Size distribution over 300s')
xlabel('Total sites in a burst')
ylabel('Count')
nFit = n ;
nFit(~isfinite(log(n))) = [] ;
xout(~isfinite(log(n))) = [] ;
p = polyfit(log(xout(25:end)),log(nFit(25:end)),1) ;
y = exp(polyval(p,log(xout))) ;
hold on
loglog(xout,y)
str = {'p = ',num2str(p(1))};
text(max(xout)-1,max(y),str)
%% Exponential fit
% subplot(1,2,1)
% Signal = Duration;
[N,edges] = histcounts(Signal,100);
Y = (edges(1:end-1)+edges(2:end))/2;
x = Y(15:end)';
y = N/sum(N);
y = y(15:end)';
f = fit(x,y,'exp1')
plot(f,Y,N/sum(N),'.')
xlabel('Duration(ms)')
ylabel('Probability')
legend('data','fit exponential')
% text(-0.28,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% subplot(1,3,2)
% Signal = Size;
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% x = Y(11:end)';
% y = N/sum(N);
% y = y(11:end)';
% f = fit(x,y,'exp1')
% plot(f,Y,N/sum(N),'.')
% xlabel('Size(a.u.)')
% ylabel('Probability')
% legend('data','fit exponential')
% subplot(1,2,2)
% Signal = Interval;
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% x = Y(11:end)';
% y = N/sum(N);
% y = y(11:end)';
% f = fit(x,y,'exp1')
% plot(f,Y,N/sum(N),'.')
% xlabel('Interval(ms)')
% ylabel('Probability')
% legend('data','fit exponential')
% text(-0.28,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%% Levy alpha-stable distribution fit
% subplot(1,2,2)
% Signal = Interval(Interval > 0); % Interval(Interval > 0); % [-Interval,Interval]; % Duration; % [-Duration,Duration]
% Signal = [-Signal,Signal];
[N,edges] = histcounts(Signal,100);
Y = (edges(1:end-1)+edges(2:end))/2;
semilogy(Y,N/sum(N),'.')
% set(gca, 'YScale', 'log')
pd = fitdist(Signal','Stable')
y = pdf(pd,Y);
hold on
semilogy(Y,y/sum(y),'m-.','LineWidth',2);
% pd = fitdist(Dis,'normal')
% y = pdf(pd,Y);
% hold on
% semilogy(Y,y/sum(y),'r--','LineWidth',2);
% legend('data','fit levy \alpha stable','Gaussian')
% xlim([0 7e2])
% ylim([6e-6 1])
xlabel('Interval(ms)')
ylabel('Probability')

% text(-0.28,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% subplot(1,3,2)
% Signal = Size;
% Signal = [-Signal,Signal];
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% semilogy(Y,N/sum(N),'.')
% set(gca, 'YScale', 'log')
% pd = fitdist(Signal','Stable')
% y = pdf(pd,Y);
% hold on
% semilogy(Y,y/sum(y),'r--','LineWidth',2);
% legend('data','fit levy \alpha stable')
% xlabel('Size(a.u.)')
% ylabel('Probability')
% subplot(1,2,2)
% Signal = Interval;
% [N,edges] = histcounts(Signal,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% semilogy(Y,N/sum(N),'.')
% set(gca, 'YScale', 'log')
% pd = fitdist(Signal','Stable')
% y = pdf(pd,Y);
% hold on
% semilogy(Y,y/sum(y),'r--','LineWidth',2);
% legend('data','fit levy \alpha stable')
% xlabel('Interval(ms)')
% ylabel('Probability')
% text(-0.28,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%% Other Fit
Signal = Duration;
[N,edges] = histcounts(Signal,100);
Y = (edges(1:end-1)+edges(2:end))/2;
z = N/sum(N);
f = fittype(@(coe,lamda,Y) 1./(1+coe*(Y).^(lamda/2)),'independent','Y','dependent','z');
[Sfit,~] = fit(Y(:),z(:),f,'StartPoint',[1 2],'Lower',[-100 0],'Upper',[1000 10]);
y = 1./(1+Sfit.coe*(Y).^(Sfit.lamda/2));
figure
plot(Y,z,'*')
hold on
loglog(Y,y,'r.-')
str = {'lamda = ',num2str(Sfit.lamda)};
xlabel('Interval(ms)')
ylabel('Probability')
text(mean(Y),mean(y),str)
% end