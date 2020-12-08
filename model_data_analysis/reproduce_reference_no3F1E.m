% LFP amplitude predict interval
No = 1;
LFPBroad = R.LFP.LFP_broad(No,:);
LFPGamma = R.LFP.LFP_gamma(No,:);
t11 = (2e4 + 1 + 44*400):(2e4 + 400 + 44*400);
t12 = (2e4 + 1 + 19*400 + 20):(2e4 + 400 + 19*400 + 20);
t21 = (4e3 + 1 + 44*80):(4e3 + 80 + 44*80);
t22 = (4e3 + 1 + 19*80 + 4):(4e3 + 80 + 19*80 + 4);
y = vec2mat(R.num_spikes{1},5);%firing rate
y = sum(y,2);
y = y/3969*2e3;

figure(1)
subplot(2,1,1)
plot([1:400],LFPBroad(t11),'k',[1:400],LFPBroad(t12),'k:',[1:400],LFPGamma(t11),'r',[1:400],LFPGamma(t12),'r:');
set(gca,'XTickLabel',[0:5:40])
xlabel('Time(ms)')
ylabel('LFP(a.u.)')
subplot(2,1,2)
[c1,~,~] = fit([1:80]',y(t21),'smoothingspline');
[c2,~,~] = fit([1:80]',y(t22),'smoothingspline');
plot(c1,'g',[1:80],y(t21),'b')
hold on;
plot(c2,'g:');
legend('off')
hold on;
plot([45 45],[-2 3.6],'k:',[62 62],[-3 7.6],'k:')
annotation('arrow',[0.45 0.56],[0.18 0.18],'LineStyle',':','HeadLength',5,'HeadWidth',5)
annotation('arrow',[0.45 0.73],[0.16 0.16],'HeadLength',5,'HeadWidth',5)
set(gca,'XTickLabel',[0:5:40])
xlabel('Time(ms)')
ylabel('Population Frequency(Hz)')