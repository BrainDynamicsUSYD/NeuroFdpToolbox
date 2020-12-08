function Visulization3DSpikesBurst(R,S,B)
% working on 10*10 electrodes
LFP = R.LFP.LFP{1};
s = size(LFP);
bin = 10;
fs = 1e4/bin;
LFP_grid = flip(reshape(LFP(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
LFP_grid = LFP_grid - mean(LFP_grid,3); % demean
timev = [];
for i = 1:length(B.WCentroids)
    timev = [timev B.WCentroids{i}(:,1)'];
end
center =B.WCentroids;
Y = prctile(S.SpikesGrid(:),95);
j = 1;
finish = center{j}(end,1);
vidObj = VideoWriter('SpikesBurst30ms1kHz95pr.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
fig = figure;
n0 = sum(sum(S.SpikesGrid(:,:,24300)));
t0 = 8200;
for t = [8200:8260 24300:24340 26680:26720 31370:31410] % 1:length(LFP_grid)
    subplot(2,2,1) % LFP + spikes
    imagesc(LFP_grid(:,:,t))
    SpikesGrid = S.SpikesGrid(:,:,t);
    ts = sprintf('time = %8.1f ms',t);
    hold on;
    [row,col] = find(SpikesGrid);
    plot(row,col,'k.')
    title(ts);
    h1 = subplot(2,2,2);
    n = sum(SpikesGrid(:));
    plot([t0 t],[n0 n],'-o')
    n0 = n;
    t0 = t;
    hold on
    if ~mod(t,50)
        cla(h1);
    end
    title('Num of Spikes');
    subplot(2,2,3)
    GGrid = zeros(sqrt(s(1)));
    GGrid(SpikesGrid >= Y) = 1;
    currentBurst = SpikesGrid.*GGrid;
    if ~ismember(t,timev) % This is for minimum length limit for bursts %
        currentBurst = zeros(sqrt(s(1)));
    end
    imagesc(currentBurst)
    title('Burst');
    h2 = subplot(2,2,4);
    while t > finish
        j = j + 1;
        finish = center{j}(end,1);
    end
    if t < center{j}(1,1)
        cla(h2)
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
end