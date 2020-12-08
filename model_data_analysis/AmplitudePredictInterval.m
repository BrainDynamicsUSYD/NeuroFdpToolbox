% LFP amplitude predict interval
No = 1;
LFPBroad = R.LFP.LFP_broad(No,:);

% R.num_spikes      2*1cell  [1*1e5 double]
%                            [1*1e5 double]
% The collection of spike numbers in two populations along the simulation time.
SpikeE = R.num_spikes{1};
window_size1 = 50; % 5 ms
window_size2 = 10; % 1 ms
SpikeE = tsmovavg(SpikeE,'s',window_size1,2);
SpikeE = tsmovavg(SpikeE,'s',window_size2,2);
SpikeI = R.num_spikes{2};
SpikeI = tsmovavg(SpikeI,'s',window_size1,2);
SpikeI = tsmovavg(SpikeI,'s',window_size2,2);
[p1,l1] = findpeaks(SpikeE); % Find local maxima: p--peak value, l--location
[p2,l2] = findpeaks(SpikeI);
[p3,l3] = findpeaks(LFPBroad);
Peak1 = p1(1:end-1);
Peak2 = p2(1:end-1);
Peak3 = p3(1:end-1);
Interval1 = l1(2:end) - l1(1:end-1);
Interval2 = l2(2:end) - l2(1:end-1);
Interval3 = l3(2:end) - l3(1:end-1);
figure
for i = 1 % :3
    subplot(11,2,[8*i-6,8*i-4,8*i-2])
    switch i
        case 1
            dat = [Peak1',Interval1'];
        case 2
            dat = [Peak2',Interval2'];
        case 3
            dat = [Peak3',Interval3'];
    end
    n = hist3(dat,[50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
    h = pcolor(xb,yb,n1);
    set(h, 'EdgeColor', 'none');
    h.ZData = ones(size(n1)) * -max(max(n));
    colormap(hot)
    oldcmap = colormap(gray);
    colormap( flipud(oldcmap) );
    colorbar
    switch i
        case 1
            [r,p] = corrcoef(Peak1',Interval1');
            text(-0.46,1.22,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        case 2
            [r,p] = corrcoef(Peak2',Interval2');
        case 3
            [r,p] = corrcoef(Peak3',Interval3');
    end
    if p(2) > 10^(-4)
        pvalue = ['p = ',num2str(p(2),'%.2f')];
    else
        pvalue = ['p < 10^{-4}'];
    end
    %     title(['r = ',num2str(r(2),'%.2f'),pvalue])
    tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
    ylabel('Interval(ms)')
    text(0.7,0.2,tp,'Units', 'Normalized','FontSize',8)
    switch i
        case 1
            xlabel('Excitatory spike rate(Hz)')
        case 2
            xlabel('Inhibitory spike rate(Hz)')
        case 3
            xlabel('LFP amplitude(a.u.)')
    end
end
% plot(Peak2,Interval2,'.')
% figure
% plot(Peak3,Interval3,'.')
for i = 1%:1e4
    t = 1e4 + i*400 + 1:1e4 + (i+1)*400;
    subplot(11,2,[13:2:21])
    plot(0.1*(t-1.04e4),LFPBroad(t),'k')
    ylabel('LFP(a.u.)')
    subplot(11,2,[1:2:9])
    plot(0.1*(t-1.04e4),SpikeI(t),0.1*(t-1.04e4),SpikeE(t))
    legend('inhibitory','excitatory')
    xlabel('time(ms)')
    ylabel('Spike rate(Hz)')
    text(-0.28,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    %     next = input(sprintf('\t Next time-frame?'));
    %     close all;
end