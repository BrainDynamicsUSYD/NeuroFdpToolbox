% function Visulization3DLFPBurst(R,B)
% % working on 20*20/40*40 LFP
% LFP = R.LFP.LFP{1};
% s = size(LFP);
% bin = 10;
% fs = 1e4/bin;
% LFP_grid = flip(reshape(LFP(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
% clear LFP
% % LFP_grid = LFP_grid - mean(LFP_grid,3); % demean
% %% Butterworth filter
% order = 4; % 4th order
% lowFreq = 30; % gamma band (default values for this function are 150-250 Hz)
% hiFreq = 80;
% Wn = [lowFreq  hiFreq]/(fs/2);
% [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
% % LFP_grid = LFP_grid(:,:,1*fs+1:end); % discard transient dynamics
% for m = 1:sqrt(s(1))
%     for j= 1:sqrt(s(1))
%         gamma_tmp = filter(b,a,LFP_grid(m,j,:));
%         gamma_tmp = gamma_tmp(:);
% %         gamma_tmp = gamma_tmp(1*fs+1:end-1*fs); % cut the head & tail after filtering
%         hil_tmp = abs(hilbert( gamma_tmp));
%         gamma_amp_grid(m,j,:) = hil_tmp; % (0.2*fs+1:end-0.2*fs); % cut the head & tail after Hilbert transform
%     end
% end
% Y = prctile(gamma_amp_grid(:),47);
% clear R
% %% baseline + 2SD
% [R] = GetBurst2(R);
% gamma_amp_grid = flip(reshape(R.LFP.LFP_gamma_hilbert_abs(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
% IsBurst = flip(reshape(R.LFP.GammaBurstEvent.is_burst(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
% clear R
%%
m1 = minmax(LFP_grid(:)');
m2 = minmax(gamma_amp_grid(:)');
% LFP_grid = LFP_grid(:,:,1.2*fs+1:end-1.2*fs);
timev = [];
for i = 1:length(B.WCentroids)
    timev = [timev B.WCentroids{i}(:,1)'];
end
center =B.WCentroids;
j = 1;
finish = center{j}(end,1);
% vidObj = VideoWriter('LFPGammaBurst30ms1kHz95pr.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 6; % number of frames to display per second
% open(vidObj);
fig = figure;
for t = 450:length(LFP_grid)
    subplot(2,2,1)
    imagesc(LFP_grid(:,:,t),m1) % pay attention to flip
    ts = sprintf('time = %8.1f ms  LFP', t);
    title(ts);
    subplot(2,2,2)
    A = gamma_amp_grid(:,:,t);
    imagesc(A,m2) % pay attention to flip
    title('Gamma Amp');
    subplot(2,2,3)
    GGrid = IsBurst(:,:,t);
%     GGrid = zeros(sqrt(s(1)));
%     GGrid(A >= Y) = 1;
    currentBurst = A.*GGrid;
    if ~ismember(t,timev) % This is for minimum length limit for bursts % t+2.2*fs
        currentBurst = zeros(sqrt(s(1)));
    end
    imagesc(currentBurst) % pay attention to flip
    title('Burst');
    h = subplot(2,2,4);
    while t > finish % t+2.2*fs
        j = j + 1;
        finish = center{j}(end,1);
    end
    if t < center{j}(1,1) % t+2.2*fs
        cla(h)
    else
        Trajectory = center{j}(:,[2:3 1]);
        
%         m = find(Trajectory(:,3)==t);
%         plot(Trajectory(m,2),Trajectory(m,1),'o')
%         set(h, 'Ydir', 'reverse')
%         set(h, 'YAxisLocation', 'Right')
%         xlim([0 40])
%         ylim([0 40])

%         [MSD,tau] = get_MSD_PBC(Trajectory(1:end,:)/40*600);
%         p = polyfit(log(tau(1:10)),log(MSD(1:10)),1);
%         loglog(tau,MSD,'.-')
%         y = exp(polyval(p,log(tau))) ;
%         hold on
%         loglog(tau,y)
%         title('Mean Square Distance versus \tau within a burst')
%         xlabel('\tau (ms)')
%         ylabel('Mean Square Distance (electrode)')
%         str = {'p = ',num2str(p(1))};
%         text(max(tau),max(y),str)
    end
%     F = getframe(fig);
%     writeVideo(vidObj,F.cdata);
    pause(0.05);
end

% load('0011-LocalSpikesGrid.mat','SpikesGrid');
% gamma_amp_grid = SpikesGrid;

% close(gcf);
% close(vidObj);

%% 63*63
% load('0011-201805261712-16526_in_1527318856057_0_neurosamp.mat','gamma_power_grid');
load('3DBurstLFP0011minTime30SR1000P95.mat','WCentroids');
bin = 10;
timev = [];
for i = 1:length(WCentroids)
    timev = [timev WCentroids{i}(:,1)'];
end
center = WCentroids;
j = 1;
finish = center{j}(end,1);
gamma_amp_grid = gamma_power_grid(:,:,1:bin:end);
s = size(gamma_amp_grid);
Y = prctile(gamma_amp_grid(:),95);
vidObj = VideoWriter('LFPGammaBurst30ms1kHz95pr.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
fig = figure;
for t = 6680:7800 % 1:length(gamma_amp_grid)
    subplot(1,3,1)
    A = gamma_amp_grid(:,:,t);
    imagesc(A)    
    ts = sprintf('time = %8.1f ms gamma amp', t);
    title(ts);
    subplot(1,3,2)
    GGrid = zeros(s(1));
    GGrid(A >= Y) = 1;
    currentBurst = A.*GGrid;
    if ~ismember(t,timev) % This is for minimum length limit for bursts % t+2.2*fs
        currentBurst = zeros(s(1));
    end
    imagesc(currentBurst)
    title('Burst');
    h = subplot(1,3,3);
    while t > finish % t+2.2*fs
        j = j + 1;
        finish = center{j}(end,1);
    end
    if t < center{j}(1,1) % t+2.2*fs
        cla(h)
    else
        Trajectory = center{j}(:,[2:3 1]);
        [MSD,tau] = get_MSD(Trajectory(1:end,:));
        p = polyfit(log(tau(1:10)),log(MSD(1:10)),1);
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
    F = getframe(fig);
    writeVideo(vidObj,F.cdata);
    pause(0.05);
end
close(gcf);
close(vidObj);
% end