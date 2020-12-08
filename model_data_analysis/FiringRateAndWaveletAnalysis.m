% function [t,ind] = FiringRateAndWaveletAnalysis(R)
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% dir_strut2 = dir('*0_neurosamp.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
H = [];
L = [];
dt = 1;
fs = 1e3; % 1e4;
lowFreq = 30;
hiFreq = 80; % Hz
Wn = [lowFreq hiFreq]/(fs/2);
order = 4; % 4th order
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
gaus_width = 12.5; % ms
[ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );
for id_out = 29 % :num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    %     S = load(files2{id_out});
    % bin = 10; % time window length
    % FR = mean(vec2mat(R.Analysis.Hz_t{1},bin),2);
    % FR = R.Analysis.Hz_t{1};
    % FR = R.num_spikes{1}/3969*1e4;
    FR = vec2mat(R.num_spikes{1},10);
    FR = sum(FR,2)'/3969*1e3;
%     FR = R.LFP.LFP{1}(3,:);
    % FR = R.grid.num_spikes_win/3969*200; % Hz
    %     I_exc = S.I_AMPA + S.I_ext; % mean(vec2mat((S.I_AMPA(1,:) + S.I_ext(1,:)),10),2);
    %     I_inh = S.I_GABA;
    %     I_tot = I_exc + I_inh;
    %     subplot(1,3,1)
    %     histogram(I_exc(:),'Normalization','probability')
    %     xlabel('I_{exc}(nA)')
    %     ylabel('probability')
    %     subplot(1,3,2)
    %     histogram(I_inh(:),'Normalization','probability')
    %     xlabel('I_{inh}(nA)')
    %     ylabel('probability')
    %     subplot(1,3,3)
    %     histogram(I_tot(:),'Normalization','probability')
    %     xlabel('I_{tot}(nA)')
    %     ylabel('probability')
    %     I_exc = mean(vec2mat(S.I_exc,10),2);
    
    
    %     all = length(FR);
%     fc = centfrq('cmor1.5-1');
%     FRcertain = filter(b,a,FR);
    %     I_exc_gamma = filter(b,a,I_exc);
    
    % baseline + 2s.d. bursts
%     hil = abs(hilbert(FRcertain));
%     gamma_hilbert_abs = conv(hil, Kernel,'same');
%     GammaBurstEvent = GetBurstFR(gamma_hilbert_abs);
    %     H = [H GammaBurstEvent.burst_du_steps{1}];
    %     L = [L GammaBurstEvent.flat_du_steps{1}];
    
    % Wavelet 95th bursts
%     scalerange = fc./([lowFreq hiFreq]*(1/fs));
%     scales = scalerange(end):0.5:scalerange(1);
%     pseudoFreq = scal2frq(scales,'cmor1.5-1',1/fs);
%     coeffs_tmp = abs(cwt(FR,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
%     CData = transpose(coeffs_tmp);
%     Y = prctile(CData(:),95); % 95
%     CData(CData < Y) = 0;
%     GreyImage = CData;
%     CData(CData>0) = 1;
%     binaryImage = CData;
%     cwt(FR,fs,'VoicesPerOctave',30);
    [wt,f,coi] = cwt(FR,fs,'VoicesPerOctave',30);
    for j = 1:length(coi)
        ind = find(f<=coi(j));
        wt(ind,j) = NaN;
    end
    ind = find(f<0 | f >80);
    wt(ind,:) = [];
    f(ind) = [];
    CData2 = abs(wt);
    Y = prctile(CData2(:),95); % 95
    CData2(CData2 < Y) = 0;
    GreyImage2 = CData2;
    CData2(CData2>0) = 1;
    binaryImage2 = CData2;

    win = 1e4-1-2.4e3; % 2e4-1
    %     t_mid = R.grid.t_mid;
    %     SP = R.grid.bayes.bayes_factor_ln >log(100);
    %     %     t_mid = 5:10:1.2e6;
    for i = 1+2.4e3:win:length(FR)-win % 60001:win:length(FR)-win
        %         %         [wt,f,coi] = cwt(FRcertain(i:i+win),fs,'TimeBandwidth',60);
        %         %         ind1 = find(f >= lowFreq & f <= hiFreq);
        %         %         ind2 = find(coi<30);
        %         %         f = f(ind1);
        %         %         wt = wt(ind1,:);
        t = i:i+win; % [i:i+win]; % *R.dt;
        %         t = round(t_mid(i:i+win)*R.dt);
        %         %         t = t(ind2);
        %         %         fr = FRcertain(i:i+win);
        %         fr = FRcertain(t);
        %         fr = fr(ind2);
        %     coeffs_tmp = abs(cwt(FRSWR(i:i+win),scales,'cmor1.5-1'))'; % this had narrow width of scales
        %     CData = transpose(coeffs_tmp); % coeffs_tmp'
        figure
%         subplot(3,1,1)
%         imagesc(abs(wt(:,t)))
%         plot(t,FRcertain(t),'b',t,hil(t),'g',t,gamma_hilbert_abs(t),'r')
%         hold on;
%         plot(t,10*GammaBurstEvent.is_burst(t),'k')
        %         plot(1e-3*t,abs(fr))
        %         xlim(1e-3*[t(1) t(end)])
        %         title('Population Firing Rate')
        %         xlabel('Time(s)')
        %         ylabel('Frequency(Hz)')
%         subplot(3,1,2)
%         imagesc(binaryImage(:,t))
%         [cfs,frq] = cwt(FRcertain(t),fs); % i:i+win
        tms = (0:numel(FR(t))-1)/fs; % FRcertain
        surface(tms,f,abs(wt(:,t))) % cfs
        axis tight
        shading flat
        xlabel('Time(s)')
        ylabel('Frequency(Hz)')
%         set(gca,'yscale','log','ytick',[4 8 16 32 64 128 256])
        %         subplot(3,1,2)
        %         t = 1e4+1:2e4;
        %         plot(0.1*t,I_exc(t),t*0.1,I_inh(t))
        %         legend('Excitatory','Inhibitory')
        %         title('Synaptic Input to One Neuron')
        %         xlabel('Time(ms)')
        %         ylabel('nA')
%         subplot(3,1,3)
%         imagesc(binaryImage2(:,t))
        %         plot(t,I_exc_gamma(t))
        %         title('Gamma-band Excitatory Input to One Neuron')
        %         xlabel('Time(ms)')
        %         xlim(t([1 end]));
        %         ylabel('Population rate(Hz)')
        
        %         figure
        %         RasterPlotYL(R,1,1,[],t)
        %         figure
        %         figure
        %         cwt(I_exc(t),fs)
        %         figure
        %         cwt(I_exc_gamma(t),fs)
        
        
        %                 figure
        %                 sp = SP(i:i+win);
        %                 bar(sp*2,'FaceColor',[0.7 0.8 0.9])
        %         plot(t,sp)
        %         inBetween = [zeros(1,length(t)), fliplr(sp)];
        %         fill([t, fliplr(t)], inBetween,[0.7 0.8 0.9]);
        %                 [wt,f,coi] = cwt(I_exc(t),fs);
        %                 for j = 1:length(coi)
        %                     ind = find(f<=coi(j));
        %                     wt(ind,j) = NaN;
        %                 end
        %                 CData = abs(wt);
        %                 Y = prctile(CData(:),95); % 95
        %                 CData(CData < Y) = 0;
        %                 ind = nansum(CData);
        %                 ind(ind>0) = 1;
        %                 bar(3*ind,'FaceColor',[0.7 0.8 0.9])
        
        %                 [wt,f,coi] = cwt(FRcertain(t),fs);
        %                 for j = 1:length(coi)
        %                     ind = find(f<=coi(j));
        %                     wt(ind,j) = NaN;
        %                 end
        %                 CData = abs(wt);
        %                 Y = prctile(CData(:),95); % 95
        %                 CData(CData < Y) = 0;
        %                 ind = nansum(CData);
        %                 %     ind = find(ind);
        %                 %     uimagesc(t,f(end:-1:1),CData)
        %
        %                 ind(ind>0) = 1;
        %                 hold on;
        %                 bar(ind,'FaceColor',[0.5 0.4 0.6])
        
        %                 [wt,f,coi] = cwt(I_exc_gamma(t),fs);
        %                 for j = 1:length(coi)
        %                     ind = find(f<=coi(j));
        %                     wt(ind,j) = NaN;
        %                 end
        %                 CData = abs(wt);
        %                 Y = prctile(CData(:),95); % 95
        %                 CData(CData < Y) = 0;
        %                 ind = nansum(CData);
        %                 ind(ind>0) = 1;
        %                 hold on;
        %                 bar(ind,'FaceColor',[0.3 0.9 0.5])
        %                 legend('Spikes Pattern','Population Fr')
        %                 legend('Excitatory Input to One Neuron','Gamma-band Population Fr','Gamma-band Excitatory Input to One Neuron')
        %         xlabel('Time(ms)')
        %         plot(t,ind)
        %         inBetween = [zeros(1,length(t)), fliplr(ind)];
        %         fill([t, fliplr(t)], inBetween,[0.5 0.4 0.6]);
        %         [~,high_du1,low_du1,high_start,~] = seq_postprocess(ind,1); % ms
        %         H = [H high_du1];
        %         L = [L low_du1];
        %
        %         uimagesc(t,f(end:-1:1),abs(wt(end:-1:1,ind2)))
        %         ylim(f([end,1]))
        %         uimagesc(R.dt*[i:i+win],pseudoFreq(end:-1:1),CData(end:-1:1,:))
        
        %         colorbar('off');
        %         ylim([lowFreq hiFreq]);
        %         xlabel('Time(ms)')
        %         ylabel('Frequency(Hz)');
        %         set(gca,'YDir','normal')
        next = input('\t Next figure?');
        close all
        %         delete(gcf);
    end
end