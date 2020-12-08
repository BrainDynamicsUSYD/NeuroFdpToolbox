figure_width = 17.8; % cm
figure_hight = 12; % cm
figure('NumberTitle','off','name', 'GammaPaperFig2', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!
% follow Fig.1
V_I = -80; % mV, reversal potential
V_E = 0;
LFPind = [1 13 13 14 14 11 12 12]; % new
%     CorrR1 = zeros(1,8);
%     CorrR2 = zeros(1,8);
tempAmp = [];
tempICI = [];
EPSC = (V_I - V_E)*(S.I_AMPA(:) + S.I_ext(:))./(S.V(:) - V_E);
IPSC = abs((V_E - V_I)*S.I_GABA(:)./(S.V(:) - V_I));
oneLFP = R.LFP.LFP{1}(LFPind,:);
dat1 = [oneLFP(:),EPSC];
dat2 = [oneLFP(:),IPSC];
n = hist3(dat1,[50 50]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat1(:,1)),max(dat1(:,1)),size(n,1)+1);
yb = linspace(min(dat1(:,2)),max(dat1(:,2)),size(n,1)+1);
subplot(2,3,1)
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
[r,p] = corrcoef(oneLFP(:),EPSC);
if p(2) > 10^(-4)
    pvalue = ['p = ',num2str(p(2),'%.2f')];
else
    pvalue = ['p < 10^{-4}'];
end
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([0 100 100 0 0],[0 0 6 6 0],'k')
xlabel('LFP(a.u)','fontsize',10)
ylabel('|EPSC(nA)|','fontsize',10)
xlim([0 100])
ylim([0 6])
text(-0.28,1,'A','Units', 'Normalized','FontSize',12)

subplot(2,3,2)
n = hist3(dat2,[50 50]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat2(:,1)),max(dat2(:,1)),size(n,1)+1);
yb = linspace(min(dat2(:,2)),max(dat2(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
[r,p] = corrcoef(oneLFP(:),IPSC);
if p(2) > 10^(-4)
    pvalue = ['p = ',num2str(p(2),'%.2f')];
else
    pvalue = ['p < 10^{-4}'];
end
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([0 100 100 0 0],[0 0 20 20 0],'k')
xlabel('LFP(a.u)','fontsize',10)
ylabel('|IPSC(nA)|','fontsize',10)
xlim([0 100])
ylim([0 20])
% text(-0.28,1,'B','Units', 'Normalized','FontSize',12)

ICI = [];
Amp = [];
for i = 1:no
    tempAmp = [];
    [p1,l1] = findpeaks(R.LFP.LFP{1}(i,:));
    [p2,l2] = findpeaks(-R.LFP.LFP{1}(i,:));
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
    tempAmp = tempAmp(1:length(tempICI));
    ICI = [ICI tempICI];
    Amp = [Amp tempAmp];
end
ind = find(ICI >= 5 & ICI <= 40);
ICI = ICI(ind);
Amp = Amp(ind);
subplot(2,3,3)
dat = [Amp',ICI'];
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
[r,p] = corrcoef(Amp',ICI');
if p(2) > 10^(-4)
    pvalue = ['p = ',num2str(p(2),'%.2f')];
else
    pvalue = ['p < 10^{-4}'];
end
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([xb(1) xb(end) xb(end) xb(1) xb(1)],[yb(1) yb(1) yb(end) yb(end) yb(1)],'k')
xlabel('Amplitude(a.u.)','fontsize',10)
ylabel('Interval(ms)','fontsize',10)
% text(-0.28,1,'C','Units', 'Normalized','FontSize',12)

ICI = [];
Spikes = [];
for i = 1:no
    [~,l1] = findpeaks(R.LFP.LFP{1}(i,:));
    tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
    tempSpikes = zeros(1,length(l1)-1);
    for j = 2:length(l1)-1
        tempSpikes(j-1) = sum(R.num_spikes{1}(l1(j):l1(j+1)));
    end
    ICI = [ICI tempICI];
    Spikes = [Spikes tempSpikes];
end
ind = find(ICI >= 5 & ICI <= 40);
ICI = ICI(ind);
Spikes = Spikes(ind);
FR = Spikes./(ICI*1e-3)/R.N(1);
ind = find(FR < 100);
ICI = ICI(ind);
FR = FR(ind);
dat = [ICI',FR'];
n = hist3(dat,[50 50]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
subplot(2,3,4)
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
[r,p] = corrcoef(ICI',FR');
if p(2) > 10^(-4)
    pvalue = ['p = ',num2str(p(2),'%.2f')];
else
    pvalue = ['p < 10^{-4}'];
end
tp = {['r = ',num2str(r(2),'%.2f')]};% ,pvalue};
text(0.5,0.8,tp,'Units', 'Normalized','FontSize',8)
hold on;
plot([xb(1) xb(end) xb(end) xb(1) xb(1)],[yb(1) yb(1) yb(end) yb(end) yb(1)],'k')
ylabel('Firing Rate(Hz)','fontsize',10)
xlabel('Duration(ms)','fontsize',10)
text(-0.28,1,'B','Units', 'Normalized','FontSize',12)

dir_strut = dir('*PPC_LFPNo1.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
PPCmat = [];
for i = 1:num_files
    load(files{i},'PPCvalue');
    PPCmat = [PPCmat;PPCvalue];
end
subplot(2,3,[5 6])
shadedErrorBar(frequency(1):frequency(2),mean(PPCmat),std(PPCmat),'lineprops','b');
xlim(frequency)
xlabel('Frequency(Hz)','fontsize',10)
ylabel('PPC value','fontsize',10)
text(-0.28,1,'C','Units', 'Normalized','FontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc GammaPaperFig2 % this is the trick!!