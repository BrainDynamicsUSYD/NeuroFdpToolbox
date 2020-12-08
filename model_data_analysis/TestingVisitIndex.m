function peak = TestingVisitIndex
%% VisitIndex
Interval = 10:5:200; % ms
ratio = 0.7:0.05:1.3;
peak = zeros(1,length(Interval));
for i = 1:length(Interval)
    I = VisitIndex(Interval(i));
    [~,ind] = max(I);
    peak(i) = ratio(ind);
end
plot(Interval,peak)
xlabel('Interval(ms)')
ylabel('Peak IE ratio')
saveas(gcf,'IntervalVsIEratio.eps')
disp('Done.')
%% VisitIndex2
% plot(0.7:0.05:1.3,cellfun(@mean,Duration),'o-')
Duration = nanmean(reshape(cellfun(@nanmean,Duration),[13,10])'); 
Displacement = nanmean(reshape(cellfun(@nanmean,Displacement),[13,10])');
%
subplot(1,3,1)
plot(0.7:0.05:1.3,Duration,'o-') % 0.7:0.05:1.3 % 0.81:0.01:0.99
xlabel('IE ratio')
ylabel('Period(ms)')
subplot(1,3,2)
plot(0.7:0.05:1.3,Displacement,'o-')
xlabel('IE ratio')
ylabel('Displacement(a.u.)')
title('12 sites')
subplot(1,3,3)
plot(0.7:0.05:1.3,Displacement./Duration,'o-')
xlabel('IE ratio')
ylabel('Speed(a.u.)')
%% plot LFP pattern duration VS interval
dir_strut = dir('3DBurst30*minTime0SR1000.mat'); % 3DBurst30*minTime0SR1000 % *RYG
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
ratio = 0.7:0.05:1.3;
for i = 1:13
    interval = [];
    duration = [];
    for j = 1:10
        R = load(files{13*(j-1)+i});
        interval = [interval R.centInterval(2:end)];
        duration = [duration R.Duration(1:end-1)];
    end
    subplot(3,5,i)
    plot(duration,interval,'.')
    xlabel('duration(ms)')
    ylabel('interval(ms)')
    title(sprintf('IE ratio=%0.2f',ratio(i)))
end
end