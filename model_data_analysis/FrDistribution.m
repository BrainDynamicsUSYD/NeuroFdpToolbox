% E firing rate distribution
% dir_strut = dir('*RYG.mat');
% num_files = length(dir_strut);
% files = cell(1,num_files);
% for id_out = 1:num_files
%     files{id_out} = dir_strut(id_out).name;
% end
edges = [0.2:0.26:8];
E = (edges(1:end-1)+edges(2:end))/2;
NN = [];
MeanStd = [];
% for i = 15 % 1:num_files
%     fprintf('Loading RYG.mat file %s...\n', files{i});
%     R = load(files{i});
    %     disp(min(R.Analysis.rate{1}))
    %     disp(max(R.Analysis.rate{1}))
    rate = R.Analysis.rate{1};
    rate = rate(rate > 0);
    [N,~] = histcounts(rate,edges);
    parmhat = lognfit(rate);
    NN = [NN;N];
    MeanStd = [MeanStd;parmhat];
% end
MS = nanmean(MeanStd,1);
plot(E,nanmean(NN,1),'o')
% err = std(NN);
% subplot(2,2,1)
% errorbar(E,nanmean(NN),err,'o');
set(gca, 'XScale', 'log');
hold on;
x = 0.1:0.1:10;
y = lognpdf(x,MS(1),MS(2));
semilogx(x,1100*y)
xlim([0.1 10])
text(-0.28,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
xlabel('Firing Rate(Hz)')
ylabel('Number of Neurons')