function [R] = getGammaAmp(R) % (stdin)
% I/EPSC balance + multi gamma features
% getting LFP Gamma ampiltude between peak to next bottom

R = GetICI(R);
LFP_gamma = R.LFP.LFP_gamma;
No = 3;
dt = R.dt;
dir_strut = dir('*0_neurosamp.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
Amp = [];

V_I = -80; % mV, reversal potential
V_E = 0;
% S = load(files{1});
% EPSC = (V_I - V_E)*(S.I_AMPA(1,:) + S.I_ext(1,:))./(S.V(1,:) - V_E);
% IPSC = (V_E - V_I)*S.I_GABA(1,:)./(S.V(1,:) - V_I);
t1 = 7e4 + 1;
t2 = 7e4 + 4e3;
t = t1:t2;
LFP_broad = R.LFP.LFP_broad(No,:);
[no,~] = size(R.LFP.LFP_gamma); % R.LFP.LFP_broad; R.pop_stats.V_mean{1}

% set(gcf,'color','w');
for i = 1:no
    [p1,l1] = findpeaks(LFP_gamma(i,:));
    [p2,l2] = findpeaks(-LFP_gamma(i,:));
    l = min(length(l1),length(l2));
    m = 1;
    for j = 1:l
        while l1(j) >= l2(m) && m <= l
            m = m + 1;
        end
        if m > l
            break
        end
        Amp = [Amp p1(j) + p2(m)];
        m = m + 1;
    end
end

R.LFP.GammaAmp = Amp;
figure('color','w')
% subplot(2,2,1)
% plot(t,IPSC(t),'b',t,4*EPSC(t),'r')
% set(gca,'XTickLabel',[0:200:400])
% ylabel('PSC(nA)')
% xlabel('Time(ms)')
% legend('IPSC','4*EPSC')
% text(-0.28,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(2,2,3)
max_lag = 50; % ms 50
[ac,lags] = autocorr(LFP_broad, round(max_lag/dt) );
plot(lags*dt,ac,'k',-lags*dt,ac,'k')
ylim([-0.2 1]) % -0.1
%     set(gca,'XTickLabel',[-50 0 50])
xlabel('Time(ms)')
ylabel('Autocorrelation')
text(-0.28,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(2,2,2)
histogram(Amp,50)
xlim([0 500]) % 100
xlabel('LFP(a.u.)')
ylabel('count')
title('Amplitude Distribution','FontWeight','normal')
text(-0.28,1.12,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(2,2,4)
histogram(R.LFP.GammaICI)
xlim([0 50])
xlabel('Interval(ms)')
ylabel('count')
title('ICI Distribution','FontWeight','normal')
text(-0.28,1.12,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
end