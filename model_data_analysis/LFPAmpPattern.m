% Collect LFP amplitude pattern dynamics on 2D scale
% only consider one pattern
% centroids + average width + time point

% visualize instantaneous 2D LFP signal
tall = length(ripple_power_grid);
mm = minmax(ripple_power_grid(:)');
Y = prctile(ripple_power_grid(:),96);
mm(1) = Y;
for i = 1:tall
    imagesc(flipud(ripple_power_grid(:,:,i)'),mm)
    pause(0.01)
end

%% moving steps
% load('0006-201804201449-05982_in_1524200018506_0_neurosamp.mat','gamma_power_grid');
bin = 1;
power_grid = ripple_power_grid(:,:,1:bin:end);
Y = prctile(power_grid(:),96);
s = size(power_grid);
hw = 31;
[Lattice,~] = lattice_nD(2,hw);
SWR2D = cell(1,2);
for t = 1:s(3)
    A = power_grid(:,:,t);
    %     A = flipud(A');
    GGrid = zeros(s(1));
    GGrid(A >= Y) = 1;
    if sum(GGrid(:)) > 0
        [~,ind] = max(A(:));
        % find hotpot centroid
        [I_row, I_col] = ind2sub(s([1 2]),ind);
        A_c = circshift(A, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
        GGrid = circshift(GGrid, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
        S = regionprops(GGrid,A_c,{'Centroid','WeightedCentroid'});
        centroids = cat(1, S.WeightedCentroid);
        pos = centroids-[round(s(1)/2)-I_row round(s(2)/2)-I_col]-32;
        %         % directly select peak as center
        %         pos = Lattice(ind,:);
        SWR2D{2} = [SWR2D{2};pos];
        SWR2D{1} = [SWR2D{1} t];
    end
end
dis = [];
for i = 1:length(SWR2D{1})
    x = SWR2D{2}(i,1);
    y = SWR2D{2}(i,2);
    if i == 1
        x0 = x;
        y0 = y;
        continue
    end
    if y0 > -11 || y0 < -21 || y > -11 || y < -21
        continue
    end
    x = SWR2D{2}(i,1);
    y = SWR2D{2}(i,2);
    d = Distance_xy(x0,y0,x,y,63);
    dis = [dis d];
    x0 = x;
    y0 = y;
end
%% Duration
sigBinary = zeros(s);
for t = 1:s(3)
    A = power_grid(:,:,t);
    GGrid = zeros(s(1:2));
    GGrid(A >= Y) = 1;
    sigBinary(:,:,t) = GGrid;
end
CC = bwconncomp(sigBinary,6);
per = [1 1 0];
CC = CC2periodic(CC,per);
Area = regionprops(CC,'Area') ;

%%
minBurstTime = 300/bin; % 30ms
count = 1;
for iBurst = 1:size(CC.PixelIdxList,1)
    currentIdx = CC.PixelIdxList{iBurst} ;
    burstTimeEnd = floor((currentIdx(end)-1)/(s(1)*s(2))) +1 ;
    burstTimeStart = floor((currentIdx(1)-1)/(s(1)*s(2))) +1 ;
    Duration2(iBurst) = burstTimeEnd-burstTimeStart+1 ;        
    if Duration2(iBurst) < minBurstTime
        continue
    end
    Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration
    Start(count) = burstTimeStart;
    End(count) = burstTimeEnd;
    count = count + 1;
end
%% MSD
slopeAll = [];
msdSuper = [];
for iBurst = 1:length(Duration)
    [~,si] = find(SWR2D{1}==Start(iBurst));
    [~,ei] = find(SWR2D{1}==End(iBurst));
        Trajectory = [SWR2D{2}(si:ei,:) SWR2D{1}(si:ei)'];
        Trajectory(:,1:2) = Trajectory(:,1:2)/63*600; %% manually modify here according to electrodes!!! %%
        [MSD,tau] = get_MSD_PBC(Trajectory);
        tempMSD = [MSD,tau] ;
        p_temp = [];
        norErr = [] ;
        maxStart = 20 ;               % minimum tau used to fit MSD
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
            slopeAll = [slopeAll pBest(1)];
            %             if pBest(1) > 1.4 && pBest(1) < 1.5  % iBurst == 9 && i == 1
            figure
            loglog(tau,MSD,'.-')
            hold on
            loglog(tau,y) % *2)
            title('Mean Square Distance versus \tau') % within a burst')
            xlabel('\tau (ms)')
            ylabel('Mean Square Distance (um)')
            str = {'p = ',num2str(pBest(1))};
            text(max(tau),max(y),str)
            next = input('\t Next figure?');
            delete(gcf);
            %             end
        else
            continue
        end
        % keep the superdiffusive one to do average
        if pBest(1)>1.1
            msdSuper = [msdSuper;tempMSD] ;
        end
end
%%
eta = 0.7; % 0.67
deltaTime = [2 4 8 16 32]; 
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;    
        for iBurst = 1:length(Duration)
            [~,si] = find(SWR2D{1}==Start(iBurst));
    [~,ei] = find(SWR2D{1}==End(iBurst));
            posCenter = SWR2D{2}(si:ei,:)/63*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
            DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            Displace = [Displace;DisplaceTemp] ;
        end
    %     rmsd = sqrt(mean(Displace));
    Displace = sqrt(Displace);
    rmsd = deltaT^eta;
    Displace = Displace/rmsd;
    [n,x] = hist(Displace,80) ; % 40
    loglog(x,n/sum(n),'o')
    hold on;
    if deltaT == 16
        pd = fitdist(Displace,'normal'); % normal
        y = pdf(pd,x) ;
        x1 = x;
    end
    if deltaT == 4
        pd = fitdist(Displace,'Stable');
        yy = pdf(pd,x);
        x2 = x;
    end
    if deltaT == 32        
        loglog(x1,y/sum(y),'r--','LineWidth',2); % 0.1*
        hold on ;
        loglog(x2,yy/sum(yy),'m-.','LineWidth',2);        
    end
    %     legend(['t = ', num2str(deltaT),' ms'],'Gaussian')
    
    %     pd = fitdist(Displace,'HalfNormal')
    %     y = pdf(pd,x) ;
    %     hold on ; loglog(x,10^3.5*y,'r--');
    %         legend('t = 30 ms','Gaussian')
    %     xlabel('r(t)')
    %     ylabel('P(r(t))')
    %     next = input('\t Next figure?');
    %     delete(gcf);
end
% xlim([0.9 1e3])
% legend('t = 2 ms','t = 4 ms','t = 8 ms','t = 16 ms','t = 32 ms','t = 40 ms')
legend('t = 2 ms','t = 4 ms','t = 8 ms','t = 16 ms','t = 32 ms','Gaussian','Stable')
xlabel('displacement/\eta^{0.8}')
ylabel('Probability')
title('Displacement Distribution')
%% Collecting Process
C = []; % centroids n*2
W = []; % average width
TP = []; % time point
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
    S = regionprops(CC,'Centroid');
    centroids = cat(1, S.Centroid);
    if count == 2
        centroids = centroids(validRegionIdx,:);
        centroids = centroids - [round(s(1)/2)-I_row  round(s(2)/2)-I_col];
        centroids = centroids - 32;
        C = [C;centroids];
        B = regionprops(CC,'BoundingBox'); % Smallest rectangle containing the region
        Boundary = cat(1, B.BoundingBox);
        W = [W mean(Boundary([3 4]))];
        TP = [TP t];
    end
end
W = 600/63*W; % um
TP = 0.1*TP; % ms
%% Visualizing Instantaneous LFP Amplitude
vidObj = VideoWriter('InstanGammaLFPAmpOver95.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30; % number of frames to display per second
open(vidObj);
for t = 200:s(3)
    A = gamma_power_grid(:,:,t);
    GGrid = zeros(s([1 2]));
    GGrid(A >= Y) = 1;
    imagesc(A.*GGrid)
    colorbar
    ts = sprintf('time = 2s + %8.1f ms', t*0.1);
    title(ts);
    pause(0.05);
    writeVideo(vidObj, getframe(gca));
    if t > 0.3e4
        break;
    end
end
close(gcf);
close(vidObj);
%% Data Analysis
% pattern size and interval
Size = []; % pattern size (each event)
Interval = [];
inds = 1;
inde = 1;
while inde < length(TP)
    if TP(inde + 1) - TP(inde) > 1
        if inde-inds >= 300
            Size = [Size max(W(inds:inde))];
            Interval = [Interval TP(inde + 1) - TP(inde)];
        end
        inds = inde + 1;
    end
    inde = inde + 1;
end