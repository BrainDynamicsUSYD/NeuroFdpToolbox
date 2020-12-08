function SSCoherence
% coherence on frequency
dir_strut = dir('*_0_neurosamp.mat');% ('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
load(files{i},'I_AMPA','I_GABA','V');
fs = 1e4; % Hz
f = [0:250]; % Hz
Cxy = [];
V_I = -80; % mV, reversal potential
V_E = 0;
EPSC = (V_I - V_E)*(I_AMPA(1,:) + I_ext(1,:))./(V(1,:) - V_E);
IPSC = (V_E - V_I)*I_GABA(1,:)./(V(1,:) - V_I);
clear I_AMPA I_GABA V
for i = 1 % :num_files
    % start form .mat files
%     fprintf('Loading RYG.mat file %s...', files{i});
%     R = load(files{i});
    disp('done.\n');
    S1 = R.LFP.LFP_broad(1,:); % gamma_hilbert_abs(1,:);
    S2 = EPSC;
    %     S2 = R.LFP.LFP_broad(2,:); %{1}(1,:);
    %     % Butterworth filter
    %     order = 4; % 4th order
    %     lowFreq = 4;
    %     hiFreq = 10;
    %     Wn = [lowFreq hiFreq]/(fs/2);
    %     [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
    %     S2 = filter(b,a,S2);
    %     S1 = angle(hilbert(S1));
    %     S2 = angle(hilbert(S2));
    %     [~,I] = sort(R.Analysis.rate{1});
    %     vS = nchoosek(I(end-49:end),2); % choose all pair combinations from the vector
%     S1 = double(full(R.spike_hist{1}(randperm(3969,1),:)));
%     for j = 1:size(R.LFP.LFP{1},1) % length(vS)
%         %         S1 = double(full(R.spike_hist{1}(vS(j,1),:)));
%         S2 = R.LFP.LFP{1}(j,:);
        [cxy,~] = mscohere(S1,S2,[],[],f,fs);
%         plot(f,cxy)
%         xlabel('Frequency(Hz)')
%         ylabel('Spike-LFP Cohenrece')
%         next = input('\t Next figure?');
%         delete(gcf);
        Cxy = [Cxy;cxy];
%     end
end
window = 10; % 3 points~1 Hz
Coh = mean(Cxy,1);
Coh = tsmovavg(Coh,'s',window,2);
figure
plot(f,Coh)
xlabel('Frequency(Hz)')
ylabel('LFP-EPSC Cohenrece')
% %% V with LFP in/out bursts
% % start form .mat files
% fprintf('Loading RYG.mat file %s...', files{1});
% R = load(files{1});
% [R] = GetBurst(R);
% disp('done.\n');
% figure
% for mode = 1:2
%     Cxy = [];
%     for i = 1:length(R.LFP.GammaBurstEvent.burst_du_steps)
%         if mode == 1
%             l = length(R.LFP.GammaBurstEvent.burst_du_steps{i});
%             start = R.LFP.GammaBurstEvent.burst_start_steps{i};
%             duration = R.LFP.GammaBurstEvent.burst_du_steps{i};
%         else
%             l = length(R.LFP.GammaBurstEvent.flat_du_steps{i});
%             start = R.LFP.GammaBurstEvent.burst_start_steps{i}+R.LFP.GammaBurstEvent.burst_du_steps{i};
%             start = start(1:end-1);
%             duration = R.LFP.GammaBurstEvent.flat_du_steps{i};
%         end
%         for j = 1:l
%             Ind = start(j):(start(j)+duration(j)-1);
%             V = R.pop_stats.V_mean{1}(Ind);
%             if length(V) < 9
%                 continue
%             end
%             LFP = R.LFP.LFP_broad(i,Ind);
%             [cxy,~] = mscohere(V,LFP,[],[],f,fs); % need at least 8 elememts in V/LFP
%             Cxy = [Cxy;cxy];
%         end
%     end
%     window = 1; % 3 points~1 Hz
%     Coh = mean(Cxy,1);
%     Coh = tsmovavg(Coh,'s',window,2);
%     plot(f,Coh)
%     hold on
% end
% legend('in-burst','out-burst')
% xlabel('Frequency(Hz)')
% ylabel('V-LFP Cohenrece')
end