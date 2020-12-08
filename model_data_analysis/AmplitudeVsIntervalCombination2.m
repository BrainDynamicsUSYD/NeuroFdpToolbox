function AmplitudeVsIntervalCombination2(R,S)
% LFP amplitude predict PREVIOUS interval
% NEED load RYG.mat and 0_neurosamp.mat
LFP_gamma = R.LFP.LFP_gamma;
t = 1e4:1e5;
V_E = 0;
V_I = -80; % mV
g_E = -nanmean((S.I_AMPA(:,t) + S.I_ext(:,t))./(S.V(:,t) - V_E));% uS
g_I = -nanmean(S.I_GABA(:,t)./(S.V(:,t) - V_I));
[no,~] = size(R.LFP.LFP_gamma); % R.LFP.LFP_broad; R.pop_stats.V_mean{1}
ICI = [];
Amp = [];

% set(gcf,'color','w');
for i = 1:no
    tempAmp = [];
    [p1,l1] = findpeaks(LFP_gamma(i,:));
    [p2,l2] = findpeaks(-LFP_gamma(i,:));
    l = min(length(l1),length(l2));
    m = 1; % index for l2
    for j = 1:l
        while m <= length(l2) && l1(j) >= l2(m)
            m = m + 1;
        end
        if m > length(l2)
            break
        end
        tempAmp = [tempAmp p1(j) + p2(m)];
        m = m + 1;
    end
    tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
    if length(tempAmp) < length(tempICI)
        disp('Miss Match!');
    end
    tempAmp = tempAmp(2:end);
    tempICI = tempICI(1:length(tempAmp));
    ICI = [ICI tempICI];
    Amp = [Amp tempAmp];
end
figure
for i = [1 3]
    subplot(2,2,i)
    switch i
        case 3
            dat = [Amp',ICI'];
        case 1
            dat = [g_E',g_I'];
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
        case 3
            [r,p] = corrcoef(Amp',ICI');
        case 1
            [r,p] = corrcoef(g_E',g_I');
    end
    if p(2) > 10^(-4)
        pvalue = ['p = ',num2str(p(2),'%.2f')];
    else
        pvalue = ['p < 10^{-4}'];
    end
    %     title(['r = ',num2str(r(2),'%.2f'),pvalue])
    tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
    text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
    switch i
        case 3
            text(-0.28,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
            xlabel('amplitude(a.u.)')
            ylabel('Interval(ms)')
        case 1
            text(-0.28,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
            xlabel('g_E(nS)')
            ylabel('g_I(nS)')
    end
end
end