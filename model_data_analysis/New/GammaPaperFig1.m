% Gamma Pattern Paper Figure 1
figure_width = 17.8; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'GammaPaperFig1', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

R = load('0001-201804161614-20473_in_1523859588189_out_RYG.mat');
S = load('0001-201804161614-20473_in_1523859588189_0_neurosamp.mat');

[R] = GetBurst(R);
dt = R.dt;
step_tot = R.step_tot;
seg_size = 1e4; % 2s
frequency = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;
pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies 80 ~ 30 Hz
i = 1;
seg = 18;
seg_ind = get_seg(step_tot, seg_size, seg);
t = (seg_ind-1)*dt*1e-3; % second
fs = 1/(dt*1e-3); % sampling frequency (Hz)
subplot(3,1,1);
x_tmp = R.LFP.LFP_gamma(i,seg_ind);
coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
CData = transpose(coeffs_tmp); % coeffs_tmp'
plot(t,R.LFP.LFP_broad(i,seg_ind),'g')% 'color',[0.8 0.8 0.8])
hold on;
plot(t,x_tmp,'b') % ,'LineWidth',1)
xlabel('Time(s)','fontsize',10);
legend('broadband','gamma-band') % legend({'broadband','gamma-band'},'fontsize',10)
text(-0.1,1,'A','Units', 'Normalized','FontSize',12) %,'FontWeight','bold')
ylabel('LFP(a.u.)','fontsize',10) % 14*LFP ('LFP(uV)')
subplot(3,1,2)
uimagesc(t,pseudoFreq(end:-1:1),CData(end:-1:1,:));
% colorbar('off');
xlim([seg_size/1e4*(seg-1) seg_size/1e4*seg]);
ylim(frequency);
xlabel('Time(s)','fontsize',10);
ylabel('Frequency(Hz)','fontsize',10);
set(gca,'YDir','normal')
set(gca,'Clipping','Off')
start = R.LFP.GammaBurstEvent.burst_start_steps{i};
ending = start-1 + R.LFP.GammaBurstEvent.burst_du_steps{i};
start = start(start>=(seg-1)*seg_size & start<=seg*seg_size);
ending = ending(ending>=(seg-1)*seg_size & ending<=seg*seg_size);
ls = length(start);
le = length(ending);
if ls == le
    l = ls;
elseif ls > le
    l = le;
    start  = start(1:l);
else
    l = ls;
    ending = ending(1:l);
end
hold on;
plot([dt*1e-3*start;dt*1e-3*start],[30*ones(1,l);140*ones(1,l)],'r--');
hold on;
plot([dt*1e-3*ending;dt*1e-3*ending],[30*ones(1,l);140*ones(1,l)],'r--');

[no,steps] = size(R.LFP.LFP_broad);
PEi = [];
IND = logical(R.LFP.GammaBurstEvent.is_burst);
for j = 1:no
    HilbertOne = hilbert(R.LFP.LFP_gamma(j,:));
    pei = repelem(angle(HilbertOne(IND(j,:))),R.num_spikes{1}(IND(j,:)));
    PEi = [PEi pei];
end
subplot(3,5,11)
polarhistogram(PEi)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
rticklabels({})
text(-0.28,1,'B','Units', 'Normalized','FontSize',12)

subplot(3,5,[12 13])
t = 5e3+(1:1e4);
I_exc = S.I_AMPA + S.I_ext;
I_inh = S.I_GABA;
plot(t*dt,I_exc(1,t),t*dt,I_inh(1,t))
legend('Excitatory','Inhibitory')
xlabel('Time(ms)','fontsize',10)
ylabel('Current(nA)','fontsize',10)
text(-0.1,1,'C','Units', 'Normalized','FontSize',12)

h = subplot(3,5,[14 15]);
p = get(h,'position');
[N,edges] = histcounts(I_exc(1,:),100);
Y = (edges(1:end-1)+edges(2:end))/2;
pd = fitdist(I_exc(1,:)','Stable')
y = pdf(pd,Y);
plot(Y,N/sum(N),'LineWidth',2)
hold on
plot(Y,y/sum(y),'m-.','LineWidth',2);
xlabel('I_{exc}(nA)','fontsize',10)
ylabel('PDF','fontsize',10)
text(-0.1,1,'D','Units', 'Normalized','FontSize',12)
ax=axes('Position',[0.75 0.2 0.1 0.1],'Unit','normalize',...
    'parent',1);
% ax=axes('Position',[1.1 1.1 .4 .4].*p,'Unit','normalize',...
%     'parent',1);
box on
loglog(ax,Y,N/sum(N),'LineWidth',2)
hold on
loglog(ax,Y,y/sum(y),'m-.','LineWidth',2);
% xlim(minmax(I_exc1(:)'))
xlim([0.5 max(I_exc(5,:)')])
set(gca,'TickLength',[0.05, 0.01])
set(gca,'XTickLabel',[],'YTickLabel',[])

% set(gca,'FontSize', 6);  % 6 points for x-axis tickmark labels
% xlabel('30 point label', 'fontsize', 30 ); % this must be after the above line!

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc GammaPaperFig1 % this is the trick!!