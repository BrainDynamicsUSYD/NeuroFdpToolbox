% analyze a series of RYG.mat
%% RYG
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%% config_data
dir_strut2 = dir('*0_neurosamp.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end
%% synaptic plasticity
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
% IndC = [513 2497];
for id_out = 19:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    %     post_dist = lattice_nD_find_dist(Lattice,hw,IndC(1)); % IndC(2)
    %     [~,IndP] = sort(post_dist);
    %     LP = R.ExplVar.local_population;
    %     A1 = IndP(1:LP)';
    %     post_dist = lattice_nD_find_dist(Lattice,hw,IndC(2)); % IndC(2)
    %     [~,IndP] = sort(post_dist);
    %     A2 = IndP(1:LP)';
    ux = (R.syn_sample.u).*(R.syn_sample.x);
    %     load(files2{id_out-1},'all');
    %     R.all = all;
    % show u & x in different zones
    figure
    subplot(5,1,1)
    plot(0.1*(1:R.step_tot),R.syn_sample.u(1,:),0.1*(1:R.step_tot),R.syn_sample.x(1,:));legend('u','x')
%     plot(0.1*(1:R.step_tot),mean(R.syn_sample.u(1:10,:)),0.1*(1:R.step_tot),mean(R.syn_sample.x(1:10,:)));legend('u','x')
    subplot(5,1,2)
    plot(0.1*(1:R.step_tot),R.syn_sample.u(2,:),0.1*(1:R.step_tot),R.syn_sample.x(2,:));legend('u','x')
% plot(0.1*(1:R.step_tot),mean(R.syn_sample.u(11:20,:)),0.1*(1:R.step_tot),mean(R.syn_sample.x(11:20,:)));legend('u','x')
    subplot(5,1,3)
    plot(0.1*(1:R.step_tot),R.syn_sample.u(3,:),0.1*(1:R.step_tot),R.syn_sample.x(3,:));legend('u','x')
%     plot(0.1*(1:R.step_tot),mean(R.syn_sample.u(21:30,:)),0.1*(1:R.step_tot),mean(R.syn_sample.x(21:30,:)));legend('u','x')
    subplot(5,1,4)
    plot(0.1*(1:R.step_tot),R.syn_sample.u(4,:),0.1*(1:R.step_tot),R.syn_sample.x(4,:));legend('u','x')
% plot(0.1*(1:R.step_tot),mean(R.syn_sample.u(31:40,:)),0.1*(1:R.step_tot),mean(R.syn_sample.x(31:40,:)));legend('u','x')
    subplot(5,1,5)
    plot(0.1*(1:R.step_tot),R.syn_sample.u(end,:),0.1*(1:R.step_tot),R.syn_sample.x(end,:));legend('u','x')
%     plot(0.1*(1:R.step_tot),mean(R.syn_sample.u(41:50,:)),0.1*(1:R.step_tot),mean(R.syn_sample.x(41:50,:)));legend('u','x')
    xlabel('Time(ms)')
    %     A1S = sum(R.spike_hist{1}(A1,:));
    %     A2S = sum(R.spike_hist{1}(A2,:));
    %     A3S = sum(R.spike_hist{1}(1:LP,:));
    %     figure
    %     subplot(3,1,1)
    %     plot(1:R.step_tot/4,sum(reshape(A1S,4,[])))
    %     subplot(3,1,2)
    %     plot(1:R.step_tot/4,sum(reshape(A2S,4,[])))
    %     subplot(3,1,3)
    %     plot(1:R.step_tot/4,sum(reshape(A3S,4,[])))
    %     disp(sum(A1S(1e4:2e4)))
    %     disp(sum(A2S(1e4:2e4)))
    %     disp(sum(A3S(1e4:2e4)))
    %     disp(sum(A1S(5.4e4:6.4e4)))
    %     disp(sum(A2S(5.4e4:6.4e4)))
    %     disp(sum(A3S(5.4e4:6.4e4)))
    % show u*x of different zones
    %     a = R.ExplVar.coefficient_K_weight*mean(R.syn_sample.u(1:8,:)).*mean(R.syn_sample.x(1:8,:));
    %     b = R.ExplVar.coefficient_K_weight*mean(R.syn_sample.u(11:20,:)).*mean(R.syn_sample.x(11:20,:));
    %     c = mean(R.syn_sample.u(21:30,:)).*mean(R.syn_sample.x(21:30,:));
    %     figure
    %     plot(0.1*(1:1e5),a); % ,0.1*(1:1e5),b,0.1*(1:1e5),c);legend('1','2','out')
    % show raster plot in different zones
    %     figure
    %     subplot(3,1,1)
    %     raster_plotCertain(R,1,1)
    %     subplot(3,1,2)
    %     raster_plotCertain(R,1,2)
    %     subplot(3,1,3)
    %     raster_plotCertain(R,1,3)
    next = input('\t Next figure?');
    close all
end
%% observe sample neuron
t0 = 5.36e4;
n0 = 0;
% for id_out = 5:num_files2
%     % start from .ygout  files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files2);
%     fprintf('\t File name: %s\n', files2{id_out});
%     load(files2{id_out});
for i = 1:10
    figure
    subplot(5,1,1)
    plot(0:5000,I_AMPA(n0+i,t0:t0+5000))
    %     plot(0:5000,I_AMPA(i,t0:t0+5000),0:5000,I_AMPA(n0+i,t0:t0+5000),0:5000,I_AMPA(2*n0+i,t0:t0+5000))
    title('I_{AMPA}')
    subplot(5,1,2)
    plot(0:5000,I_ext(n0+i,t0:t0+5000))
    %     plot(0:5000,I_ext(i,t0:t0+5000),0:5000,I_ext(n0+i,t0:t0+5000),0:5000,I_ext(2*n0+i,t0:t0+5000))
    title('I_{ext}')
    subplot(5,1,3)
    plot(0:5000,I_GABA(n0+i,t0:t0+5000))
    %     plot(0:5000,I_GABA(i,t0:t0+5000),0:5000,I_GABA(n0+i,t0:t0+5000),0:5000,I_GABA(2*n0+i,t0:t0+5000))
    title('I_{GABA}')
    subplot(5,1,4)
    plot(0:5000,V(n0+i,t0:t0+5000))
    %     plot(0:5000,V(i,t0:t0+5000),0:5000,V(n0+i,t0:t0+5000),0:5000,V(2*n0+i,t0:t0+5000))
    title('V')
    subplot(5,1,5)
    I_EI = I_AMPA+I_ext+I_GABA;
    plot(0:5000,I_EI(n0+i,t0:t0+5000))
    %     plot(0:5000,I_EI(i,t0:t0+5000),0:5000,I_EI(n0+i,t0:t0+5000),0:5000,I_EI(2*n0+i,t0:t0+5000))
    title('EIfight')
    next = input('\t Next figure?');
    close all
end
% end

%% Burst drift frequency
f = cell(1,num_files);
for id_out = 1:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    [Dfre,Duration,Start,End] = BurstDrifting(R);
    f{id_out} = [Dfre{:}];
end
save('DriftFrequency.mat','f')
%% firing rate
% % [g_mu,~] = CollectVectorYG('ExplVar','ExplVar.g_mu');
% % g_mu = g_mu([27:31]);
% [g_balance,~] = CollectVectorYG('ExplVar','ExplVar.g_balance');
% g_balance = g_balance([27:32]);
% [fr,~] = CollectVectorYG('Analysis','mean(Analysis.rate{1})');
% fr = fr([27:32]);
% plot(g_balance,fr,'o-')
% % xlabel('EE synaptic conductance(uS)')
% xlabel('g_{balance}')
% ylabel('Excitatory Neuron Firing Rate(Hz)')
%% plot raster plot
% figure_visibility = 'on';
% for r_num = [27 28 31 32]
%     R = load(files{r_num});
%     step_tot = R.reduced.step_tot;
%     comments = R.comments;
%     Num_pop = 1;
%     seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
%     seg_num = ceil(step_tot/seg_size);
%     for seg = 1:seg_num
%         for pop_ind = 1:Num_pop
%             % Plot raster plot
%             figure
%             raster_plot(R, pop_ind, seg, [] );
%             xlabel('Time(s)')
%         end
%     end
% end
%% plot 2D spikes snapshot
% hw = 31; % R.ExplVar.hw;
% [Lattice, ~] = lattice_nD(2, hw);
% t = 500;
% for r_num = 32 % [27 28 31 32]
%     R = load(files{r_num});
%     R = get_grid_firing_centre(R);
%     t_mid = R.grid.t_mid;
%     ind_a_vec = R.grid.ind_ab(1,:);
%     ind_b_vec = R.grid.ind_ab(2,:);
%     x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
%     y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
%     while t < 1e4
%         ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
%         figure('Name','Vis','color','w','NumberTitle','off');
%         plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
%         axis equal;
%         box on;
%         set(gca,'xtick',[],'ytick',[]);
%         xlim([-hw hw]);
%         ylim([-hw hw]);
%         ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);
%         title(ts);
%         next = input('\t Next figure?');
%         close all
%         t = t + 10;
%     end
% end