% reproduce reference NO.3 Fig1D/Fig2D
% load RYG.mat sample0.mat
[no,~] = size(R.LFP.LFP{1});
Vr_I = -80; % mV
Vr_E = 0; % mV
for j = 1:no
    % LFP = nanmean(R.LFP.LFP{1},1); % R.LFP.LFP{1}; R.pop_stats.V_mean{1}
%     LFP = R.LFP.LFP{1}(j,:);
    
    EPSC = (Vr_I - Vr_E)*(S.I_AMPA(1,:) + S.I_ext(1,:))./(S.V(1,:) - Vr_E);
    IPSC = (Vr_E - Vr_I)*S.I_GABA(1,:)./(S.V(1,:) - Vr_I);
    
    for t = 5e3:2e2:R.step_tot
        plot([1:2e2],abs(IPSC(t:(t + 199))),'b',[1:2e2],5*EPSC(t:(t + 199)),'r')
        set(gca,'XTickLabel',[0:5:20])
        legend('|IPSC|','5*EPSC')
        xlabel('Time(ms)')
        ylabel('PSC(nA)')
        next = input(sprintf('\t Next time-frame?'));
        close all;
    end
end
    
%     t = (2e4 + 1):(2e4 + 120);
%     %     t1 = (2e4 + 1):(2e4 + 400);
%     %     t2 = (4e3 + 1):(4e3 + 80);
%     %     y = vec2mat(R.num_spikes{1},5);%firing rate
%     %     y = sum(y,2);
%     %     y = y/3969*2e3;
%     for i = 285:665
%         h = figure('numbertitle','off','name','check_current_stats','color','w');
%         ax(1) = subplot(2,1,1);
%         %         plot([1:400],LFP(t1 + 400*i),'k');
%         plot([1:120],LFP(t + 120*i),'k')
%         set(gca,'XTickLabel',[0:2:12])
%         ylabel('LFP(a.u.)')
%         title(['j = ',num2str(j),'   i = ',num2str(i)])
%         ax(2) = subplot(2,1,2);
%         plot([1:120],abs(IPSC(t + 120*i)),'b',[1:120],5*EPSC(t + 120*i),'r')
%         %         [c,~,~] = fit([1:80]',y(t2 + 80*i),'smoothingspline');
%         %         plot(c,'g',[1:80],y(t2 + 80*i),'.');
%         %         set(gca,'XTickLabel',[0:5:40])
%         set(gca,'XTickLabel',[0:2:12])
%         legend('|IPSC|','5*EPSC')
%         xlabel('Time(ms)')
%         % ylim([0 20])
%         ylabel('PSC(nA)')
%         % linkaxes(ax,'x');
%         next = input(sprintf('\t Next time-frame?'));
%         close all;
%     end
% end