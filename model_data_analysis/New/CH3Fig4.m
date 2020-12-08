figure_width = 13; % cm
figure_hight = 9; % cm
figure('NumberTitle','off','name', 'CH3Fig4', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

h = subplot(2,3,1);
p = get(h,'position');
Amp = zeros(size(R.LFP.LFP_ripple));
for i = 1:16
    Amp(i,:) = abs(hilbert(R.LFP.LFP_ripple(i,:)));
end
Signal = Amp(:);
histogram(Signal,40,'Normalization','probability')
xlabel('Amplitude(a.u.)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.25,1,'A','Units', 'Normalized','FontSize',12)
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
ax=axes('Position',[0.25 0.8 0.1 0.1],'Unit','normalize',...
    'parent',1);
semilogx(ax,Y,N/sum(N),'.')
parmhat = lognfit(Signal);
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(ax,Y,YY/sum(YY),'r')

h = subplot(2,3,2);
p = get(h,'position');
Signal = 0.1*[R.LFP.ripple_event.ripple_du_steps{:}];
histogram(Signal,40,'Normalization','probability')
xlabel('Duration(ms)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.25,1,'B','Units', 'Normalized','FontSize',12)
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
ax=axes('Position',[0.5 0.8 0.1 0.1],'Unit','normalize',...
    'parent',1);
semilogx(ax,Y,N/sum(N),'.')
parmhat = lognfit(Signal);
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(ax,Y,YY/sum(YY),'r')

h = subplot(2,3,4);
p = get(h,'position');
Signal = 0.1*[R.LFP.ripple_event.flat_du_steps{:}];
histogram(Signal,40,'Normalization','probability')
xlabel('Interval(ms)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.25,1,'C','Units', 'Normalized','FontSize',12)
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
ax=axes('Position',[0.25 0.33 0.1 0.1],'Unit','normalize',...
    'parent',1);
semilogx(ax,Y,N/sum(N),'.')
parmhat = lognfit(Signal);
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(ax,Y,YY/sum(YY),'r')

% LFP Spike phase lock inside bursts
PEi = [];
Spikes = [];
[no,steps] = size(R.LFP.LFP_broad);
IND = logical(R.LFP.ripple_event.is_SWR);
for j = 1:no
    ind1 = R.LFP.ripple_event.ripple_start_steps{j};
    ind2 = R.LFP.ripple_event.ripple_start_steps{j}+R.LFP.ripple_event.ripple_du_steps{j}-1;
    for i = 1:length(ind1)
        Spikes = [Spikes sum(R.num_spikes{1}(ind1(i):ind2(i)))];
    end
    HilbertOne = hilbert(R.LFP.LFP_ripple(j,:));
    pei = repelem(angle(HilbertOne(IND(j,:))),R.num_spikes{1}(IND(j,:)));
    PEi = [PEi pei];
end
h = subplot(2,3,5);
p = get(h,'position');
Signal = Spikes;
histogram(Signal,40,'Normalization','probability')
xlabel('Spikes','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.25,1,'D','Units', 'Normalized','FontSize',12)
[N,edges] = histcounts(Signal,60);
Y = (edges(1:end-1)+edges(2:end))/2;
ax=axes('Position',[0.5 0.33 0.1 0.1],'Unit','normalize',...
    'parent',1);
semilogx(ax,Y,N/sum(N),'.')
parmhat = lognfit(Signal);
hold on
YY = lognpdf(Y,parmhat(1),parmhat(2));
semilogx(ax,Y,YY/sum(YY),'r')

subplot(1,3,3)
polarhistogram(PEi)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
rticklabels({})
% title('Phase Locking in Bursts')
text(-0.25,1,'E','Units', 'Normalized','FontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig4 % this is the trick!!