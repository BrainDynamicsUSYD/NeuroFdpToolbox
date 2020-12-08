function LocalSpikesSignal
% get local firing rate signal to plot power spectrum

% config = dir('0001*_config_data.mat');
% load(config.name,'LFP_centre_x','LFP_centre_y');
hw = 31;
LFP_centre_x = linspace(-hw, hw, 21); % E16:9  E100:21 E400:41
LFP_centre_y = linspace(-hw, hw, 21);
LFP_centre_x = LFP_centre_x(2:2:20); % E16(2:2:8)  E100(2:2:20) E400(2:2:40)
LFP_centre_y = LFP_centre_y(2:2:20);
[LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
LFP_centre_x = LFP_centre_x(:);
LFP_centre_y = LFP_centre_y(:);
dir_strut = dir('*_out_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
S1 = 25; % 20 for default
S2 = 8;
n = 1e4;
f = 1000/n*(0:round(n/2)); % start from 0.1 Hz
Local1 = {};
for i = 1:length(LFP_centre_x)
    Local1{i} = find(lattice_nD_find_dist(Lattice,hw,LFP_centre_x(i),LFP_centre_y(i)) <= S1);
    %     Local1(i,:) = find(lattice_nD_find_dist(Lattice,hw,LFP_centre_x(i),LFP_centre_y(i)) <= S1);
    %     Local2(i,:) = find(lattice_nD_find_dist(Lattice,hw,LFP_centre_x(i),LFP_centre_y(i)) <= S2);
end
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    for j = 1:length(LFP_centre_x)
        L1fr(j,:) = sum(R.spike_hist{1}(Local1{j}(:),:))/length(Local1{j})*1e3;
        %         L1fr(j,:) = sum(R.spike_hist{1}(Local1(j,:),:));
        %         L2fr(j,:) = sum(R.spike_hist{1}(Local2(j,:),:));
    end
    y = reshape(full(L1fr),[100 10 R.step_tot/10]); % spikes in ms
    y = squeeze(sum(y,2));
    yy = smooth(y(1,500:1500),'rlowess');
    findchangepts(yy,'Statistic','rms','MaxNumChanges',10)
    ipt = findchangepts(yy,'Statistic','rms','MaxNumChanges',10)
%     ydft = fft(y,n,2);
%     ydft = ydft(:,1:n/2+1);
%     psdy = (1/(1e3*n)) * abs(ydft).^2;
%     psdy(:,2:end) = 2*psdy(:,2:end);
%     for j = 1:length(LFP_centre_x)
%         %         subplot(4,4,j)
%         window_size = 40; % 4 Hz
%         psdy(j,:) = tsmovavg(psdy(j,:),'s',window_size,2);
%         plot(f(50:round(end/2)),10*log10(psdy(j,50:round(end/2))))
%         %         title('Periodogram Using FFT')
%         xlabel('Frequency (Hz)')
%         ylabel('Power/Frequency (dB/Hz)')
%         break
%     end
    next = input('\t Next figure?');
    delete(gcf);
end
end