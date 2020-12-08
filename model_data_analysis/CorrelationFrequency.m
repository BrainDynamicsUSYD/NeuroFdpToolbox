function CorrelationFrequency(stdin)
% plot power spectrum based on firing rate

if nargin == 0
    dir_strut = dir('*RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for i = 1:num_files
        files{i} = dir_strut(i).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end
rate = [];
for i = [1:10] % num_files [1 12 17 18] 
    R = load(files{i});
    pop = 1;
    %%%%%%%%%DON'T CHANGE%%%%%%%%%%
    if pop == 1
        pop_type = 'E';
    else
        pop_type = 'I';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = vec2mat(R.num_spikes{pop},10); % firing rate
    y = sum(y,2);
    y = y';
    rate = [rate;y/3969*1e3]; % Hz
end
% Butterworth filter broad band
fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
order = 2; % 2th order
lowFreq_br = [0.5:10.5 10:55 51.5:192.5];
hiFreq_br = [3.5:13.5 16:61 66.5:207.5];
l1 = length(rate);
l2 = length(lowFreq_br);
H = zeros(num_files,l1);
RHOs = zeros(1,l2);
RHOn = zeros(1,l2);
for j = 1:l2
    Wn = [lowFreq_br(j) hiFreq_br(j)]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    for i = 1:num_files
        H(i,:) = hilbert(filter(b,a,rate(i,:)));
    end
    Habs = abs(H);
    Hmean1 = mean(Habs(1:5,:));
    Hmean2 = mean(Habs(6:10,:));
    Hnoise1 = Habs(1,:) - Hmean1;
    Hnoise2 = Habs(6,:) - Hmean2;
    [RHOs(j),~] = corr(Hmean1',Hmean2','type','Spearman');
    [RHOn(j),~] = corr(Hnoise1',Hnoise2','type','Spearman');
end
save('sn3.mat','RHOs','RHOn')
end