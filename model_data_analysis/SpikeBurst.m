% spike burst distribution in E pop
% calculation start after 1s; interval < 8 ms, spikes >= 3
% Burst event rate; burst index
N = R.N(1);
BurstNumber = zeros(1,N);
BurstIndex = zeros(1,N);
tic;
for i = 1:N
    s = movsum(R.spike_hist{1}(i,:),80);
    s = s(41+1e4:end-39);
    BurstNumber(i) = sum(s>=3);
    BurstIndex(i) = sum((s>=3).*s)/sum(R.spike_hist{1}(i,:));
end
toc;
BurstRate = BurstNumber/(1e-4*(R.step_tot-1e4));
[NR,edgesR] = histcounts(BurstRate,30);
ER = (edgesR(1:end-1)+edgesR(2:end))/2;
[NI,edgesI] = histcounts(BurstIndex,30);
EI = (edgesI(1:end-1)+edgesI(2:end))/2;
subplot(1,2,1)
semilogx(NR,ER)
subplot(1,2,2)
semilogx(NI,EI)