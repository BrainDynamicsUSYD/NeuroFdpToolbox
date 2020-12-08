function LFPAmpPatternAnalysis
%% 150ms:the actual min time is 15 ms
dir_strut = dir('3DBurst*minTime30SR1000.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%% Size of Bursts v.s. Time Interval between Bursts
for i = 1:num_files
    R = load(files{i});
    PS = R.patternScale(1:end-1);
    ESInterval = 0.1*R.centInterval(2:end); % end-start interval
    dat = [PS',ESInterval'];
    n = hist3(dat,[50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
    figure(1)
    subplot(2,3,i)
    h = pcolor(xb,yb,n1);
    set(h, 'EdgeColor', 'none');
    h.ZData = ones(size(n1)) * -max(max(n));
    colormap(hot)
    oldcmap = colormap(gray);
    colormap( flipud(oldcmap) );
    colorbar
    [r,p] = corrcoef(PS',ESInterval');
    if p(2) > 10^(-4)
        pvalue = ['p = ',num2str(p(2),'%.4f')];
    else
        pvalue = ['p < 10^{-4}'];
    end
    tp = {['r = ',num2str(r(2),'%.4f')],pvalue};
    text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
    xlabel('Pattern Size')
    ylabel('Interval among Bursts(ms)')
    % Distribution of Burst Size
    figure(2)
    subplot(2,3,i)
    [N,edges] = histcounts(R.patternScale,100);
    Y = (edges(1:end-1)+edges(2:end))/2;
    loglog(Y,N,'o')
    xlabel('Pattern Size')
    ylabel('Count')
    % Distribution of Propagation Step within Bursts
    if ~isfield(R,'PStep')
        PStep = [];
        for j = 1:length(R.Centroids)
            centroids = R.Centroids{j};
            for k = 1:length(centroids)-1
                PStep = [PStep Distance_xy(centroids(k,2),centroids(k,3),centroids(k+1,2),centroids(k+1,3),20)];
            end
        end
        save(files{i},'PStep','-append')
    else
        PStep = R.PStep;
    end
    figure(3)
    subplot(2,3,i)
    histogram(30*PStep,100) % um
    xlabel('Moving Distance per Step(um)')
    ylabel('Count')
    figure(4)
    subplot(2,3,i)
    [N,edges] = histcounts(30*PStep,100);
    Y = (edges(1:end-1)+edges(2:end))/2;
    loglog(Y,N,'o')
    xlabel('Moving Distance per Step(um)')
    ylabel('Count')
    % Distribution of Propagation Step within Bursts + Jumping Distance between Bursts
    AllMove = 30*[PStep R.distCent(2:end)]; % um
    figure(5)
    subplot(2,3,i)
    histogram(AllMove,100)
    xlabel('Moving Distance(um)')
    ylabel('Count')
    figure(6)
    subplot(2,3,i)
    [N,edges] = histcounts(AllMove,100);
    Y = (edges(1:end-1)+edges(2:end))/2;
    loglog(Y,N,'o')
    xlabel('Moving Distance(um)')
    ylabel('Count')
    % Distribution of Burst Jumping Distance
    figure(7)
    subplot(2,3,i)
    histogram(30*R.distCent(2:end),100)
    xlabel('Jumping Distance(um)')
    ylabel('Count')
    % Distribution of Burst Interval
    figure(8)
    subplot(2,3,i)
    histogram(ESInterval,100)
    xlabel('Burst Interval(ms)')
    ylabel('Count')
    % Size of Bursts v.s. Bursts Jumping Distance
    ESDist = 30*R.distCent(2:end); % end-start distance
    dat = [PS',ESDist'];
    n = hist3(dat,[50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
    figure(9)
    subplot(2,3,i)
    h = pcolor(xb,yb,n1);
    set(h, 'EdgeColor', 'none');
    h.ZData = ones(size(n1)) * -max(max(n));
    colormap(hot)
    oldcmap = colormap(gray);
    colormap( flipud(oldcmap) );
    colorbar
    [r,p] = corrcoef(PS',ESDist');
    if p(2) > 10^(-4)
        pvalue = ['p = ',num2str(p(2),'%.4f')];
    else
        pvalue = ['p < 10^{-4}'];
    end
    tp = {['r = ',num2str(r(2),'%.4f')],pvalue};
    text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
    xlabel('Pattern Size')
    ylabel('Jumping Distance between Bursts(um)')
end
%% MSD
for i = 1:num_files
    R = load(files{i});
    % select only 1 burst
    center = R.WCentroids ;
    numB = length(center);
    figure(10)
    subplot(2,3,i)
    for iBurst = 20 % 20 % 1:numB
        Trajectory = [] ;
        temp = zeros(length(center{iBurst}),1) ;
        TrajectoryTemp = center{iBurst}(:,[2:3 1]);
        Trajectory = [Trajectory;TrajectoryTemp] ;
        
        [MSD,tau] = get_MSD(Trajectory(1:end,:)) ;%,grid_size);
                loglog(tau,MSD,'.-')
        p = polyfit(log(tau(1:10)),log(MSD(1:10)),1) ;
                y = exp(polyval(p,log(tau))) ;
                hold on
                loglog(tau,y)
                title('Mean Square Distance versus \tau within a burst')
                xlabel('\tau (ms)')
                ylabel('Mean Square Distance (electrode)')
                str = {'p = ',num2str(p(1))};
                text(max(tau),max(y),str)
%                 text(max(tau)+5,max(y)+5,str)
        
        slope(iBurst) = p(1) ;
    end
%     histogram(slope)
%     disp(mean(slope))
%     xlabel('Slope')
%     ylabel('Count')
end
%% Visualizing Instantaneous LFP Amplitude
load('3DBurst0001minTime30SR1000.mat','Centroids');
timev = [];
for i = 1:length(Centroids)
    timev = [timev Centroids{i}(:,1)'];
end
R = load('0001-201805101630-05734_in_1525933960856_out_RYG.mat','LFP');
% LFP_grid = flip(reshape(R.LFP.LFP{1},20,20,[]));
% dt = 0.1;
% fs = 1/(dt*1e-3); % sampling frequency (Hz)
LFP_grid = flip(reshape(R.LFP.LFP{1}(:,1:10:end),20,20,[]));
fs = 1e3;
% Butterworth filter
order = 4; % 4th order
lowFreq = 30; % gamma band (default values for this function are 150-250 Hz)
hiFreq = 80;
Wn = [lowFreq  hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
gamma_power_grid = zeros(size(LFP_grid));
for i = 1:20
    for j= 1:20
        gamma_tmp = filter(b,a,LFP_grid(i,j,:));
        gamma_tmp = gamma_tmp(:);
        hil_tmp = abs(hilbert( gamma_tmp));
        gamma_power_grid(i,j,:) = hil_tmp;
    end
end
clear R LFP_grid
Y = prctile(gamma_power_grid(:),95);
%%
vidObj = VideoWriter('th30fs1e3GammaBurstAmpOver95.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
for vt = (1:3e2)+2040
    A = gamma_power_grid(:,:,vt);
    GGrid = zeros(20);
    GGrid(A >= Y) = 1;
    currentBurst = A.*GGrid;
    if ~ismember(vt,timev)
        currentBurst = zeros(20);
    end
    imagesc(currentBurst)
    colorbar
    ts = sprintf('time = %8.1f ms', vt);
    title(ts);
    pause(0.05);
    writeVideo(vidObj, getframe(gca));
end
close(gcf);
close(vidObj);
end