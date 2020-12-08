% useful code
ax1 = subplot(2,1,1);
plot(t,LFP_broad(seg_ind),'color',[0.8 0.8 0.8])
ax2 = subplot(2,1,2);
colorbar('off');
set(ax2, 'XLim', get(ax1, 'XLim'));

set(gca,'XTickLabel',[0:200:1000])

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
for i = [1]
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
end

hinfo_inh = hdf5info('0001-201608271253-02473_in_1472266547975_1_neurosamp.h5');
I_AMPA_e = hdf5read(hinfo_exc.GroupHierarchy.Datasets(1));
up = 90*ones(size(LFP_peak));
down = 85*ones(size(LFP_peak));
plot([LFP_peak;LFP_peak],[up;down],'k','LineWidth',2)

f = 1000/steps*(0:round(steps/2) + 1);
Y = fft(LFP_gamma(No,:).^2,steps);
Pyy = Y.*conj(Y)/steps;
Pyy = Pyy(1:round(steps/2) + 2);
window_size = 200;
simple = tsmovavg(Pyy,'s',window_size,2);
plot(f(1:end),simple)
set(gca,'YScale','log');
xlim([20 100])
xlabel('Hz')
ylabel('Power spectral density')

SpikeE = R.num_spikes{1};
window_size1 = 50; % 5 ms
window_size2 = 10; % 1 ms
SpikeE = tsmovavg(SpikeE,'s',window_size1,2);
SpikeE = tsmovavg(SpikeE,'s',window_size2,2);
SpikeI = R.num_spikes{2};
SpikeI = tsmovavg(SpikeI,'s',window_size1,2);
SpikeI = tsmovavg(SpikeI,'s',window_size2,2);

rowIndex = repmat((1:10)',[1 bin]);
[~, outx] = sort(rand(10,bin),2);
[~, outy] = sort(rand(10,bin),2);
shindx = sub2ind([10,bin],rowIndex,outx);
shindy = sub2ind([10,bin],rowIndex,outy);
Xsh = X(shindx);
Ysh = Y(shindy);

plot([40 40],[0 100],'--r')

aP = diff(SegmentP);
aP(aP ~= 1) = 0;
b = diff([0 find(diff(aP)) numel(aP)]);
if aP(1) == 1
    start = 1;
else
    start = 0;
end
durationP = 2*(1 + b(start:2:end)); % ms

l = length(ts);
seq = [1]; % define the first one element
for i = 2:l-1
    left = round(ts(i)-ts(i-1));
    right = round(ts(i+1)-ts(i));
    if left == 1 && right == 1
        seq = [seq 1];
    else
        if left == 1
            seq = [seq 1 0];
            continue
        end
        if right == 1
            seq = [seq 0 1];
            continue
        end
        seq = [seq 0 1 0];
    end
end
