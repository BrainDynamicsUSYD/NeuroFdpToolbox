function BurstDuration
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
duration = [];
bin = 10;
for i = 1 % :num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    %     R = GetBurst(R);
    %     d = R.dt*[R.LFP.GammaBurstEvent.burst_du_steps{:}];
    LFP = R.LFP.LFP{1}(:,1:bin:end);
    [no,l] = size(LFP);
%     try
%         LFP_gamma = R.LFP.LFP_gamma;
%     catch
        % Butterworth filter
        order = 4; % 4th order
        lowFreq = 30; % gamma band
        hiFreq = 80;
        fs = 1e4;
        Wn = [lowFreq hiFreq]/(fs/2);
        [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
        LFP_gamma = zeros(size(LFP)); 
        for j = 1:no
            LFP_gamma(j,:) = filter(b,a,LFP(j,:));
        end
%     end
    gamma_amp_grid = zeros(size(LFP));
    for j = 1:no
        gamma_amp_grid(j,:) = abs(hilbert(LFP_gamma(j,:)));
    end
    Y = prctile(gamma_amp_grid(:),95);
    for j = 1:no
        temp = zeros(1,l);
        temp(gamma_amp_grid(j,:) >= Y) = 1;
        [~,d,~,~,~] = seq_postprocess(temp,1,1);
        duration = [duration d];
    end
    % save('duration.mat','duration')
    % subplot(2,2,4)
    histogram(0.1*duration)
    title(['mean = ',num2str(mean(0.1*duration)),'ms'])
    % text(-0.28,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    xlabel('Duration(ms)')
    ylabel('Count')
end