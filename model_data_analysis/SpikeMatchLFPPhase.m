% Spike phase match to LFP (Phase lock)
tic;
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end

%% whole phase locking of E/I spikes
PE = [];
PI = [];
for i = 1 % :num_files %  13:13:num_files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    
    [no,steps] = size(R.LFP.LFP_broad);
    try
        LFP_gamma = R.LFP.LFP_gamma;
    catch
        fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
        % Butterworth filter
        order = 4; % 4th order
        lowFreq_br = 30;
        hiFreq_br = 80;
        Wn = [lowFreq_br hiFreq_br]/(fs/2);
        [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
        for j = 1:no
            LFP_gamma(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
        end
    end
    for j = 1:no
        Hilbert(j,:) = hilbert(LFP_gamma(j,:));
    end
    Hilbertm = mean(Hilbert);
    pei = repelem(angle(Hilbertm),R.num_spikes{1});
    pi = repelem(angle(Hilbertm),R.num_spikes{2});
    PE = [PE pei];
    PI = [PI pi];
end

subplot(1,2,1)
polarhistogram(PE)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
text(-1.4,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% pax = gca;
% pax.GridColor = 'r';
subplot(1,2,2)
polarhistogram(PI)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
% text(-0.46,1.12,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% title('Excitatory neuron coupling to LFP')

%% LFP Spike lock inside/outside bursts
PEi = [];
PEo = [];
figure
for i = 1 % :num_files %  13:13:num_files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    [R] = GetBurst(R);
    [no,steps] = size(R.LFP.LFP_broad);
    IND = logical(R.LFP.GammaBurstEvent.is_burst);
    try
        LFP_gamma = R.LFP.LFP_gamma;
    catch
        fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
        % Butterworth filter
        order = 4; % 4th order
        lowFreq_br = 30;
        hiFreq_br = 80;
        Wn = [lowFreq_br hiFreq_br]/(fs/2);
        [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
        for j = 1:no
            LFP_gamma(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
        end
    end
    for j = 1:no
        HilbertOne = hilbert(LFP_gamma(j,:));
        pei = repelem(angle(HilbertOne(IND(j,:))),R.num_spikes{1}(IND(j,:)));
        peo = repelem(angle(HilbertOne(~IND(j,:))),R.num_spikes{1}(~IND(j,:)));
        PEi = [PEi pei];
        PEo = [PEo peo];
    end
end
subplot(1,2,1)
polarhistogram(PEi)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
rticklabels({})
title('Phase Locking in Bursts')
text(-0.28,1.12,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(1,2,2)
polarhistogram(PEo)
thetaticklabels({'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'})
rticklabels({})
title('Phase Locking out of Bursts')
text(-0.28,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% for i = 1:num_files
%     fprintf('Loading RYG.mat file %s...', files{i});
%     R = load(files{i});
%     disp('done.\n');
%     R = get_grid_firing_centre(R);
%     R = get_IEI(R);
%
%     LFP_gamma = R.LFP.LFP_gamma;
%     binsize = 36; % bin size = 10 degrees
%     histo1 = zeros(1,binsize);
%     [no, steps] = size(R.LFP.LFP_broad);
%     fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
%     pop = 1; % ATTENTION HERE: decide before running!
%     LFP_valley_t = []; % time steps
%     hi = [];
%     LFP_Gamma = nanmean(LFP_gamma);
%
%     %%%%%%%%%DON'T CHANGE%%%%%%%%%%
%     if pop == 1
%         pop_type = 'E';
%     else
%         pop_type = 'I';
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     for j = 2:(steps - 1)
%         if LFP_Gamma(j) < LFP_Gamma(j - 1) && LFP_Gamma(j) < LFP_Gamma(j + 1)
%             valley = j;
%             LFP_valley_t = [LFP_valley_t valley];
%         end
%     end
%     l = length(LFP_valley_t);
%     step_bin = round((LFP_valley_t(2:l) - LFP_valley_t(1:(l-1)))/binsize);
%     for j = 1:binsize
%         histo1(j) = sum(R.num_spikes{pop}((LFP_valley_t(1:(l - 1)) + (j - 1)*step_bin)));
%     end
%     histo1 = histo1/sum(histo1);
%     histo(i,:) = histo1;
%     %     figure('color','w')
% %     figure
% %     bar([histo1,histo1],1)
% %     xlim([1 72])
% %     set(gca,'XTick',[1 36 72])
% %     set(gca,'XTickLabel',[10 360 720])
% %     xlabel('Gamma Phase(^o)')
% %     ylabel(['Spike Probability(',pop_type,')'])
% end
% histo = mean(histo);
% figure
%     bar([histo,histo],1)
%     xlim([1 72])
%     set(gca,'XTick',[1 36 72])
%     set(gca,'XTickLabel',[10 360 720])
%     xlabel('Gamma Phase(^o)')
%     ylabel(['Spike Probability(',pop_type,')'])