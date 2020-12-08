function p = WMtoThetaPhase
% method for detecting retrieval: 80% within 20 ms
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
dir_strut2 = dir('*_config_data.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
load(dir_strut2.name,'StiNeu');
% pick = 29;
% load(files2{pick-11},'StiNeu')
NumP = length(StiNeu);
x = cell(1,NumP);
p = cell(1,NumP);
bin = 20; % ms
Color = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0];
seg = 14.5e3; % ms
% Butterworth filter
order = 4; % 4th order
lowFreq = 4; % gamma band
hiFreq = 7;
fs = 1e3;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
for id_out = 1:num_files
    % For PBS array job
    %     loop_num = loop_num + 1;
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    R.num_spikes{1} = sum(full(R.spike_hist{1}));
    FR = vec2mat(R.num_spikes{1},10);
    FR = sum(FR,2)'/3969*1e3;
    FR_theta = filter(b,a,FR);
    FR_theta_phase = angle(hilbert(FR_theta));
    tstart = 2251; % ms
    ystart = min(FR_theta);
    ystop = max(FR_theta);
    for i = 1:NumP
        r = movsum(full(R.reduced.spike_hist{1}(StiNeu{i},tstart:end)),bin,2); % Hz (id_out,:)
        r(r > 1) = 1;
        %     r{id_out} = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu{id_out},period)),length(StiNeu{id_out}),bin,[]),2))/bin*1e4; % Hz
        x{i} = tstart-1 + find(sum(r) > 40);
        p{i} = [p{i} FR_theta_phase(x{i})];
        tstart = tstart + 2e3;
    end
    %     plot(tstart:tstart+seg,FR_theta(tstart:tstart+seg))
    %     for i = 1:NumP
    %         xhere = x{i}(x{i}>=tstart & x{i}<=tstart+seg);
    %         l = length(xhere);
    %         hold on;
    %         plot([xhere;xhere],[ystart*ones(1,l);ystop*ones(1,l)],'color',Color(i,:)) % *ones(1,l)
    %     end
    %     xlabel('Time(ms)')    
%     next = input('\t Next figure?');
%     close all
    
    %     RasterPlotYL2(R,StiNeu)
    %     Top = max(FR_theta);
    %     plot(100*FR_theta-30+id_out*50)
    
end
figure
for i = 1:NumP
    subplot(2,4,i)
    polarhistogram(p{i})
end
end