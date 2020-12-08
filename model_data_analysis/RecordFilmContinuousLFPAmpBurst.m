% function RecordFilmContinuousLFPAmpBurst
% R = load('0016-201807271136-13419_in_1532655627546_out_RYG.mat','LFP');
load('0015-201806291103-38833_in_1530234329865_0_neurosamp.mat','LFP_grid');
R.LFP.LFP = {reshape(LFP_grid,63^2,[])};
% load('3DBurstLFP0016minTime30SR1000P95.mat','WCentroids');
load('3DBurst30015minTime30SR1000.mat','WCentroids');
bin = 10;
timev = [];
for i = 1:length(WCentroids)
    timev = [timev WCentroids{i}(:,1)'];
end
center = WCentroids;
%% For 95pr visualization
LFPr = R.LFP.LFP{1};
fw = sqrt(size(LFPr,1));
LFPr = flip(reshape(LFPr(:,1:bin:end),fw,fw,[]));
% LFPr = LFPr - mean(LFPr,3); % demean
% Butterworth filter
order = 4; % 4th order
lowFreq = 30; % gamma band (default values for this function are 150-250 Hz)
hiFreq = 80;
fs = 1e3;
Wn = [lowFreq  hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass');
% LFPCut
% LFPr = LFPr(:,:,1*fs+1:end);
climLFPr = minmax(LFPr(:)');
gamma_amp_grid = zeros(size(LFPr));
% gamma_amp_grid = zeros(fw,fw,size(LFPr,3)-2.4*fs);
for i = 1:fw
    for j = 1:fw
        gamma_tmp = filter(b,a,LFPr(i,j,:));
        gamma_tmp = gamma_tmp(:);
        % LFPCut
        %         gamma_tmp = gamma_tmp(1*fs+1:end-1*fs);
        gamma_amp_grid(i,j,:) = abs(hilbert(gamma_tmp));
        % LFPCut
        %         hil_tmp = abs(hilbert(gamma_tmp));
        %         gamma_amp_grid(i,j,:) = hil_tmp(0.2*fs+1:end-0.2*fs); % total cut:first2.2s+last1.2s
    end
end
clear R gamma_grid gamma_tmp
s = size(gamma_amp_grid);
Y = prctile(gamma_amp_grid(:),95);
clim = minmax(gamma_amp_grid(:)');
% vidObj = VideoWriter('PropagationGammaBurst1.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 6; % number of frames to display per second
% open(vidObj);
for num = 1:length(WCentroids)
    a = WCentroids{num}(1,1);
    b = WCentroids{num}(end,1);
    period = a:b;
    j = 1;
    finish = center{j}(end,1);
    fig = figure;
    set(gcf,'Position',[681 707 560 254])
    % period = 1190:420;
    temp =gamma_amp_grid(:,:,period);
    mm = minmax(temp(:)');
    for t = period % 1:length(gamma_amp_grid)
        A = gamma_amp_grid(:,:,t);
        subplot(1,2,1)
        %%% 63*63
        %     MemP = V(:,:,t);
        %     imagesc(flipud(MemP),climV)
        %     title('Mem Potential')
        %%%% 20*20
        %         LFP = LFPr(:,:,t+1.2*fs);
        %     imagesc(flipud(LFP),climLFPr)
        imagesc(flipud(A),clim)
        subplot(1,2,2)
        
        GGrid = zeros(s(1));
        GGrid(A >= Y) = 1;
        currentBurst = A.*GGrid;
        %         if ~ismember(t,timev) % This is for minimum length limit for bursts % t+2.2*fs
        %             currentBurst = zeros(s(1));
        %         end
        imagesc(flipud(currentBurst),clim)
        ts = sprintf('Burst, t = %8.1f ms',t);
        title(ts);
        %     while t > finish % t+2.2*fs
        %         j = j + 1;
        %         finish = center{j}(end,1);
        %     end
        %     h1 = subplot(2,2,3);
        %     xlim([0 40])
        %     ylim([0 40])
        %     if ismember(t,center{j}(:,1))
        %         if t == center{j}(1,1)
        %             plot(center{j}(1,3),center{j}(1,2),'k.','MarkerSize',3)
        %             temp = center{j}(1,2:3);
        %             ind = 1;
        %         else
        %             ind = ind + 1;
        %             y = [temp(1) center{j}(ind,2)];
        %             x = [temp(2) center{j}(ind,3)];
        %             hold on
        %             plot(x,y,'-k.','LineWidth',0.5,'MarkerSize',3)
        %             temp = center{j}(ind,2:3);
        %         end
        %     else
        %         cla(h1)
        %     end
        %     h2 = subplot(2,2,4);
        %     if t < center{j}(1,1) % t+2.2*fs
        %         cla(h2)
        %     else
        %         Trajectory = center{j}(:,[2:3 1]);
        %         Trajectory(:,1:2) = (Trajectory(:,1:2) - 0.5)/63*600;
        %         [MSD,tau] = get_MSD_PBC(Trajectory(1:end,:));
        %         p_temp = [];
        %         norErr = [] ;
        %         maxStart = 20 ;               % minimum tau used to fit MSD
        %         tau_max = maxStart:min(100,size(MSD,1)) ;
        %         for fitIdx = 1:length(tau_max)
        %             fitRange = (1:tau_max(fitIdx)) ;
        %             [pAll,~] = polyfit(log(tau(fitRange)),log(MSD(fitRange)),1) ;
        %             y = exp(polyval(pAll,log(tau(fitRange)))) ;
        %             errorRate = mean(abs((MSD(fitRange)-y(fitRange))./MSD(fitRange))) ;
        %             p_temp(fitIdx,:) = pAll ;
        %             norErr(fitIdx) = errorRate ;
        %         end
        %         bestIdx = find(norErr == min(norErr)) ;
        %         pBest = p_temp(bestIdx(end),:) ;
        %         y = exp(polyval(pBest,log(tau))) ;
        %         loglog(tau,MSD,'.-')
        %         hold on
        %         loglog(tau,y)
        %         title('Mean Square Distance versus \tau within a burst')
        %         xlabel('\tau (ms)')
        %         ylabel('Mean Square Distance (electrode)')
        %         str = {'p = ',num2str(pBest(1))};
        %         text(max(tau),max(y),str)
        %     end
        %             F = getframe(fig);
        %             writeVideo(vidObj,F.cdata);
        pause(0.2);
    end
    next = input('\t Next figure?');
    delete(gcf);
end
% close(gcf);
% close(vidObj);
%% Burst propagation snapshots 95pr
j = 1;
finish = center{j}(end,1);
fig = figure;
for i = 1:8
    t = 7856 + 5*i;
    subplot(2,4,i)
    A = gamma_amp_grid(:,:,t);
    GGrid = zeros(s(1));
    GGrid(A >= Y) = 1;
    currentBurst = A.*GGrid;
    if ~ismember(t,timev) % This is for minimum length limit for bursts % t+2.2*fs
        currentBurst = zeros(s(1));
    end
    imagesc(flipud(currentBurst),clim)
    ts = sprintf('t = %8.1f ms',t);
    title(ts);
    while t > finish % t+2.2*fs
        j = j + 1;
        finish = center{j}(end,1);
    end
    if i == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    end
end
%% For baseline+2SD visualization
[R] = GetBurst2(R);
sigBinary = R.LFP.GammaBurstEvent.is_burst(:,1:bin:end);
s = size(sigBinary);
sigBinary = reshape(sigBinary,sqrt(s(1)),sqrt(s(1)),[]);
LFP_gamma_hilbert_abs = reshape(R.LFP.LFP_gamma_hilbert_abs(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]);
clear R
s = size(sigBinary);
CC = bwconncomp(sigBinary,6);
per = [1 1 0];
CC = CC2periodic(CC,per);
minBurstTime = 300/bin;
clim = minmax(LFP_gamma_hilbert_abs(:)');
%% 
GT = [];
for iBurst = 1: size(CC.PixelIdxList,1)
    currentIdx = CC.PixelIdxList{iBurst} ;
    burstTimeEnd = floor((currentIdx(end)-1)/(s(1)*s(2))) +1 ;
    burstTimeStart = floor((currentIdx(1)-1)/(s(1)*s(2))) +1 ;
    Duration =  burstTimeEnd-burstTimeStart+1 ;
    if Duration < minBurstTime
        continue
    end
    burstIdxTemp = zeros(s) ;
    burstIdxTemp(currentIdx) = 1 ;
    currentBurst = sigBinary.*LFP_gamma_hilbert_abs.*burstIdxTemp ;
    set(gcf,'Position',[681 707 560 254])
    for iTime = burstTimeStart:burstTimeEnd
        a = find(currentBurst(32,32,:));
        GT = [GT;a];
        A = LFP_gamma_hilbert_abs(:,:,iTime);
        subplot(1,2,1)
        imagesc(flip(A),clim)
        subplot(1,2,2)
        Burst = currentBurst(:,:,iTime);
        imagesc(flip(Burst),clim)
        ts = sprintf('Pattern, t = %8.1f ms',iTime);
        title(ts);
        pause(0.2);
    end
    next = input('\t Next figure?');
    delete(gcf);
end
%% Pattern propagation snapshots baseline+2SD
fig = figure;
for i = 1:8
    t = 210 + 5*i;
    subplot(2,4,i)    
    Burst = currentBurst(:,:,t);
    if ~ismember(t,timev) % This is for minimum length limit for bursts % t+2.2*fs
        Burst = zeros(s(1));
    end
    imagesc(Burst,clim)
    ts = sprintf('t = %8.1f ms',t);
    title(ts);
    if i == 1
        text(-0.2,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    end
end
%% Burst trajectory
j = 6;
h1 = figure;
% h1 = subplot(3,4,[9 10]);
for t = 908:991
    if ismember(t,center{j}(:,1))
        if t == center{j}(1,1)
            plot(63-center{j}(1,2),center{j}(1,3),'k.','MarkerSize',5)
            temp = [63-WCentroids{j}(1,2),WCentroids{j}(1,3)];
            ind = 1;
        else
            ind = ind + 1;
            y = [temp(2) center{j}(ind,3)];
            x = [temp(1) 63-center{j}(ind,2)];
            hold on
            plot(x,y,'-k.','LineWidth',0.5,'MarkerSize',5)
            temp = [63-WCentroids{j}(ind,2),WCentroids{j}(ind,3)];
        end
    else
        cla(h1)
    end
    xlim([0 63])
    ylim([0 63])
end
title('Pattern Trajectory')
text(-0.1,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%% MSD
Trajectory = [center{j}(:,[2:3 1])];
Trajectory(:,1:2) = (Trajectory(:,1:2) - 0.5)/40*600; %% manually modify here according to electrdoes!!! %%
[MSD,tau] = get_MSD_PBC(Trajectory);
p_temp = [];
norErr = [] ;
maxStart = 10 ;               % minimum tau used to fit MSD
tau_max = maxStart:min(100,size(MSD,1)) ;
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
    y = exp(polyval(pBest,log(tau))) ;
    figure
    %     subplot(3,4,[11 12])
    loglog(tau,MSD,'.-')
    hold on
    loglog(tau,y)
    title('Mean Square Distance versus \tau') % within a burst')
    xlabel('\tau (ms)')
    ylabel('Mean Square Distance (um)')
    str = {'p = ',num2str(pBest(1))};
    text(max(tau),max(y),str)
    text(-0.1,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
end
% end