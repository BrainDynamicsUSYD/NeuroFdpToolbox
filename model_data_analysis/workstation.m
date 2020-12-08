% % function workstation(stdin)
% % if nargin == 0
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end

% % % % plot coherence between E/I firing rate
% % Cxy_mat = [];
% % %
% % % Cxy_mats = [];
% clear IE_ratio fr ratio ie A
[fr,~] = CollectCellYG('Analysis','mean(Analysis.rate{1})');
[ie,~] = CollectCellYG('neuron_stats','nanmean(neuron_stats.IE_ratio{1})');
ie = [ie{:}];
tefr = (te+530)./fr;
ie = vec2mat(ie,18);
ie = mean(ie);
te = vec2mat(te,18);
tefr = vec2mat(tefr,18);
error1 = std(te);
error2 = std(tefr);
te = mean(te);
tefr = mean(tefr);
[ie,I] = sort(ie);
subplot(2,1,1)
errorbar(ie,te(I),error1(I))
hold on;
plot([-0.76 -0.76],[-2 4],'--r')
xticks([-1 -0.76 -0.6:0.2:0])
xlabel('IE ratio')
ylabel('Transfer Entropy(TE)')
subplot(2,1,2)
errorbar(ie,tefr(I),error2(I))
hold on;
plot([-0.94 -0.94],[-2 4],'--r')
xticks([-1 -0.94 -0.8:0.2:0])
xlabel('IE ratio')
ylabel('TE/firing rate')
% figure
% errorbar(ie,tefr2(I),error3(I))
% for i = 1:num_files % [1 12 17 18]
    
%     % start form .mat files
%     fprintf('Loading RYG.mat file %s...', files{i});
%     R = load(files{i});
%     E = load(files2{i});
%     I = load(files3{i});
%     disp('done.\n');
%     a = nanmean(E.I_GABA,2)./(nanmean(E.I_AMPA,2) + nanmean(E.I_ext,2));
%     IE_ratio(i) = nanmean(a);
%     totE1 = E.I_AMPA + E.I_GABA + E.I_ext;
%     totI1 = I.I_AMPA + I.I_GABA + I.I_ext;
%     EE = mean(mean(E.I_AMPA,2))/mean(mean(totE1,2)); % totE2
%     IE = mean(mean(E.I_GABA,2))/mean(mean(totE1,2)); % totE2
%     EI = mean(mean(E.I_AMPA,2))/mean(mean(totI1,2)); % totE2
%     II = mean(mean(E.I_GABA,2))/mean(mean(totI1,2)); % totE2
%     A(i) = IE*EI-EE*II;
%     b = R.neuron_stats.IE_ratio{1}(R.neuron_sample.neuron_ind{1});
%     ratio(i) = nanmean(b);
%     ie(i) = nanmean(R.neuron_stats.IE_ratio{1});
%     %     R = get_IEI(R);
%     %     %
%     %     %     y1 = vec2mat(R.num_spikes{1},10); % firing rate
%     %     %     y1 = sum(y1,2);
%     %     %     y1 = y1';
%     %     %     y1 = y1/3969*1e3;
%     %     %     y1s = y1(randperm(length(y1)));
%     %     %
%     %     %     %     y2 = vec2mat(R.num_spikes{2},10); % firing rate
%     %     %     %     y2 = sum(y2,2);
%     %     %     %     y2 = y2';
%     %     %     %     y2 = y2/3969*1e3;
%     %     %     Cxy_m = [];
%     %     %     Cxy_ms = [];
%     %     %     for j = 1:16
%     %     %         y2 = R.LFP.LFP_gamma(j,10:10:end);
%     %     %         y2s = y2(randperm(length(y2)));
%     %     %         [Cxys,~] = mscohere(y1s,y2s,[],[],[],1000);
%     %     %
%     %     %         [Cxy,F] = mscohere(y1,y2,[],[],[],1000); % calculating coherence
%     %     %
%     %     %         Cxy_m = [Cxy_m Cxy];
%     %     %         Cxy_ms = [Cxy_ms Cxys];
%     %     %     end
%     %     %     n = round(length(F)/5); % Fmax=500Hz here
%     %     %     window_size = 40; % ~10 Hz
%     %     %     %     simple = tsmovavg(Cxy','s',window_size,2);
%     %     %     simple = tsmovavg(nanmean(Cxy_m,2)','s',window_size,2);
%     %     %     simples = tsmovavg(nanmean(Cxy_ms,2)','s',window_size,2);
%     %     %     Cxy_mat = [Cxy_mat;simple(1:n)];
%     %     %     Cxy_mats = [Cxy_mats;simples(1:n)];
%     %     % end
%     %     % figure
%     %     % plot(F(1:n)',mean(Cxy_mat,1),F(1:n)',mean(Cxy_mats,1));
%     %     % legend('data','shuffuled data')
%     %     % xlabel('Frequency(Hz)')
%     %     % ylabel('Coherence')
%     %     % % saveas(gcf,'Coherence_EI_fr.pdf')
%     %     % % saveas(gcf,'Coherence_Spiking_LFP.pdf')
%     %     % end
%     %
%     %
%     %     % % Plot distribution of gamma burst on frequency(first running TimeFrequencyCImage)
% % %         f = 55:(-1):5;
% % %         t = 0:100:5000;
% % %         l = length(f) - 1;
% % %         c = zeros(l);
% % %         du = R.LFP.Gamma_duration;
% % %         fr = R.LFP.Gamma_Fre;
% % %         for i = 1:l
% % %             range = du(fr >= f(i+1) & fr < f(i));
% % %             c(i,:) = histcounts(range,t);
% % %         end
% % %         imagesc(c)
% % %         colorbar
% % %         set(gca,'XtickLabel',[0.1:0.1:0.5]);
% % %         set(gca,'YtickLabel',[50:(-5):5]);
% % %         xlabel('Duration(s)')
% % %         ylabel('Frequency(Hz)')
%     %
%     %
%     %     % plot LFP amplitute-frequency
%     %     % running get_IEI in advance
%     %     peak = zeros(1,1000);
%     %     Ind = cumsum(R.LFP.Gamma_IEI_Ind + 1);
%     %     ind = [0 Ind(1:(end-1))] + 1;
%     %     amp = R.LFP.Gamma_LFP_peak;
%     %     amp(ind) = [];
%     %     fre = 1e3./R.LFP.Gamma_IEI2;
%     %     for f = 1:1000
%     %         cand = find(fre>=(f*0.1-0.05) & fre<(f*0.1+0.05));
%     %         if ~isempty(cand)
%     %             peak(f) = mean(amp(cand));
%     %         end
%     %     end
%     %     window_size = 40; % ~4 Hz
%     %     simple = tsmovavg(peak,'s',window_size,2);
%     %     Cxy_mat = [Cxy_mat;simple];
%     % end
%     % simple = tsmovavg(mean(Cxy_mat,1),'s',window_size,2);
%     % plot([1:1000]*0.1,simple)
%     % xlabel('Frequency(Hz)')
%     % ylabel('Amplitute(a.u.)')
% end
%
% % t = 1e4 + (1:1e3);
% % time = 0.1*(1:1e3);
% % for i = 1:16
% %     subplot(4,4,i)
% %     plot(time,R.LFP.LFP_gamma(i,t))

% plot firing rate-I_total function
% x = R.pop_stats.I_input_mean{2};
% y = R.Analysis.Hz_t{2};
% scatter(x,y)
% pause(1)
% end

% figure('Name','Vis','color','w','NumberTitle','off');
% % axis equal;
% % box on;
% xlim([-3 1]);
% ylim([1 3.5]);
% hold on;
% nFrames = 2000;
% vidObj = VideoWriter('FiringRate-ItotalFunction.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 7;
% open(vidObj);
% for i = 1:1e5
%     scatter(x(i),y(i))
%     ts = sprintf('time = %8.1f ms', i*0.1);
%     title(ts);
%     xlabel('Total input current(nA)')
%     ylabel('firing rate of E-pop(Hz)')
%     hold on;
%     pause(0.01)
% %     if mod(i,10) == 0
% %         writeVideo(vidObj, getframe(gca));
% %     end
% %     if  i >= 5e3
% %         break;
% %     end
% end
% close(gcf);
% close(vidObj);


% [dist] = lattice_nD_find_dist(Lattice,hw,center);
% [d,ind] = sort(dist);

