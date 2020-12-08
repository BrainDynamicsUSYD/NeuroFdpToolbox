% all kinds of color images for the collecting data
[f1,~] = CollectVectorYG('Analysis','nanmean(sqrt(Analysis.CV2_ISI{1}))');
[f2,~] = CollectVectorYG('Analysis','nanmean(sqrt(Analysis.CV2_ISI{2}))');
[f3,~] = CollectVectorYG('Analysis','nanmean(Analysis.rate{1})');
[f4,~] = CollectVectorYG('Analysis','nanmean(Analysis.rate{2})');
% % [v,loop_num] = CollectVectorYG('LFP','LFP.wavelet.peak.rp_wl_amp');
% %[v,loop_number] = CollectVectorYG('grid','nanvar(grid.jump_dist)');
% Mat = zeros(3,3);
% Hz_overall = zeros(1,6);
% CV = zeros(1,6);
% fr_E = zeros(1,6);
% fr_I = zeros(1,6);
% for sec_num = 1:6
% %     f_ripple = v(loop_num <= 5*sec_num);
% %     loop_num(loop_num <= 5*sec_num) = 1000;
% %     Mat(sec_num) = nanmean([f_ripple{:}]);
%     Hz_overall(sec_num) = nanmean(f1(sec_num + 12):27:(sec_num + 120)));%EE:(9*sec_num - 4):27:(9*sec_num + 104)
%     CV(sec_num) = nanmean(f2((sec_num + 12):27:(sec_num + 120)));%IE:(3*sec_num + 2):27:(3*sec_num + 110)
%     fr_E(sec_num) = nanmean(f3((sec_num + 12):27:(sec_num + 120)));
%     fr_I(sec_num) = nanmean(f4((sec_num + 12):27:(sec_num + 120)));
% end
subplot(2,2,1)
bar([3:11 14 17 23],f1(1:12))
% % Mat = Mat';
% imagesc(Mat)
% colorbar
% set(gca,'XTick',[1:3])
% set(gca,'XTickLabel',[3:11 14 17])
% set(gca,'YTick',[1:3])
% set(gca,'YTickLabel',[19:21])
xlabel('decay_{GABA}')
% ylabel('tau_cI')
title('CV_{ISI} of Excitatory')
% title('Firing rate of Excitatory neurons(Hz)')
subplot(2,2,2)
bar([3:11 14 17 23],f2(1:12))
% set(gca,'XTickLabel',[3:11 14 17])
xlabel('decay_{GABA}')
title('CV_{ISI} of Inhibitory')
subplot(2,2,3)
bar([3:11 14 17 23],f3(1:12))
% set(gca,'XTickLabel',[3:11 14 17])
xlabel('decay_{GABA}')
title('firing rate of Excitatory(Hz)')
subplot(2,2,4)
bar([3:11 14 17 23],f4(1:12))
% set(gca,'XTickLabel',[3:11 14 17])
xlabel('decay_{GABA}')
title('firing rate of Inhibitory(Hz)')

% saveas(gcf,'CV_ISI_firing_rate.pdf')

% % saveas(gcf,'ImageMatrix_Excitatory_CV_ISI.pdf')


% a = 7:9;
% b = 9:11;
% c = 19:21;
% d = ones(3,3);
% [A,B] = meshgrid(a,b);
% [~,C] = meshgrid(a,c);
% mesh(A,B,19*d,'FaceColor','none');
% xlabel('tau_cEE')
% ylabel('tau_cIE')
% zlabel('tau_cI')
% title('Firing rate of Inhibitory Neurons(Hz)')
% hold on;
% mesh(A,B,20*d,'FaceColor','none');
% hold on;
% mesh(A,B,21*d,'FaceColor','none');
% hold on;
% scatter3(tau_c_EE,tau_c_IE,tau_c_I,6*Hz_overall,Hz_overall,'filled')
% caxis([4,26])
% colorbar