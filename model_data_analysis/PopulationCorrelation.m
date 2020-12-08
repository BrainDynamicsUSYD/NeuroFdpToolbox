function PopulationCorrelation(StiNeu)
% Ref: "Stable population coding for working memory coexists with
% heterogeneous neural dynamics in prefrontal cortex" PNAS
bin = 25e2; % 250ms
% single file
%% no-overlap bin
steps = size(R.spike_hist{1},2)/bin;
r = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu,:)),length(StiNeu),bin,[]),2))/bin*1e4; % Hz
state1 = r(:,2e4/bin+1);
state2 = r(:,3e4/bin-1); % squeeze(sum(reshape(full(R.spike_hist{1}(Neu,29341-bin:29340)),length(Neu),bin,[]),2))/bin*1e4; % Hz
c1 = zeros(1,steps);
c2 = zeros(1,steps);
cm = zeros(steps);
for t = 1:steps
    covm = cov(state1,r(:,t));
    c1(t) = covm(2)/sqrt(var(state1)*var(r(:,t)));
    covm = cov(state2,r(:,t));
    c2(t) = covm(2)/sqrt(var(state2)*var(r(:,t)));
    for tt = 1:steps
        covm = cov(r(:,t),r(:,tt));
        cm(t,tt) = covm(2)/sqrt(var(r(:,t))*var(r(:,tt)));
    end
end
figure
plot(0.1*bin*(1:steps),c1,0.1*bin*(1:steps),c2)
xlim([2000 3000])
legend('state 1','state 2')
xlabel('Time(ms)')
ylabel('Correlation')
figure
imagesc(cm(2e4/bin-1:3e4/bin-1,2e4/bin-1:3e4/bin-1))
xticks([0.5:5.5])
yticks([0.5:5.5])
xticklabels(0.25*(7:12))
yticklabels(0.25*(7:12))
xlabel('Time(s)')
ylabel('Time(s)')


%% overlap-bin
period = 2e4+1:2.925e4;
r = movsum(full(R.spike_hist{1}(StiNeu,period)),bin,2)/bin*1e4; % Hz
state1 = r(:,0.125e4);
state2 = r(:,0.8e4); % squeeze(sum(reshape(full(R.spike_hist{1}(Neu,29341-bin:29340)),length(Neu),bin,[]),2))/bin*1e4; % Hz
steps = length(period);
c1 = zeros(1,steps);
c2 = zeros(1,steps);
% cm = zeros(steps);
for t = 1:steps
    covm = cov(state1,r(:,t));
    c1(t) = covm(2)/sqrt(var(state1)*var(r(:,t)));
    covm = cov(state2,r(:,t));
    c2(t) = covm(2)/sqrt(var(state2)*var(r(:,t)));
%     for tt = 1:steps
%         covm = cov(r(:,period(t)),r(:,period(tt)));
%         cm(t,tt) = covm(2)/sqrt(var(r(:,period(t)))*var(r(:,period(tt))));
%     end
end
figure
% c1 = smooth(c1,150);
% c2 = smooth(c2,150);
plot(0.1*(period),c1,0.1*(period),c2)
xlim([2.125e4 2.8e4])
legend('state 1','state 2')
xlabel('Time(ms)')
ylabel('Correlation')
% figure
% imagesc(cm)
% xlabel('Time(ms)')
% ylabel('Time(ms)')


%% multiple files
NumP = length(StiNeu);
r = cell(1,NumP);
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
SampleN = randperm(length(StiNeu{1}),NumP);
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    for i = 1:NumP
        % population rate
%         groupr = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu{i},:)),length(StiNeu{i}),bin,[]),2))/bin*1e4; % Hz
%         r{i} = [r{i};mean(groupr)];
        % neuron rate
        sampler = sum(reshape(full(R.spike_hist{1}(StiNeu{i}(SampleN(i)),:)),bin,[]))/bin*1e4; % Hz
        r{i} = [r{i};sampler];
    end
end
t = 0.1*bin*(1:length(R.spike_hist{1})/bin); % ms

figure
for i = 1:NumP
    meanr = mean(r{i});
    subplot(NumP,2,2*i-1)
    plot(t(1:8),meanr(1:8))
    ylim([0 10])
    xlabel('Time(ms)')
    ylabel('Rate(Hz)')
    subplot(NumP,2,2*i)
    plot(t(10:end),meanr(10:end))
    ylim([0 10])
    xlabel('Time(ms)')
    ylabel('Rate(Hz)')
end
end
