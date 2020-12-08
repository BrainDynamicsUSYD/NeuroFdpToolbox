% function FiringRate
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% dir_strut2 = dir('*_config_data.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
FR = zeros(1,130); % num_files);
for i = 1:130 % num_files
    fprintf('Processing output file No.%d out of %d...\n', i, num_files);
    fprintf('\t File name: %s\n', files{i});
    R = load(files{i});
%     load(files2{i},'StiNeu')
%     FR(i) = mean(R.Analysis.rate{1});
%     FR(i) = sum(sum(R.spike_hist{1}(StiNeu{1},2.25e4+1:end)))/50/17.75; % Hz
    FR(i) = sum(R.num_spikes{1})/3969/10; % 20s or 10s
end
FR = vec2mat(FR,13);
fr = mean(FR);
frSTD = std(FR);
Fr = 0.7:0.05:1.3;
v1 = polyfit(log10(Fr(1:7)),log10(fr(1:7)),1);
x1 = Fr(1:9);
y1 = 10^v1(2)*x1.^v1(1);
v2 = polyfit(log10(Fr(7:end)),log10(fr(7:end)),1);
x2 = Fr(5:end);
y2 = 10^v2(2)*x2.^v2(1);
errorbar(Fr,fr,frSTD,'o','MarkerSize',8,'CapSize',8,'LineWidth',1.5)
set(gca,'YScale','log')
hold on
semilogy(x1,y1,'LineWidth',2)
hold on
semilogy(x2,y2,'LineWidth',2)
xlabel('Scaled IE ratio')
ylabel('Firing Rate(Hz)')
axes('Position',[.5 .5 .35 .35])
box on
errorbar(Fr,fr,frSTD,'o','MarkerSize',8,'CapSize',8,'LineWidth',1.5)
hold on
plot(x1,y1,'LineWidth',2)
hold on
plot(x2,y2,'LineWidth',2)
-v1(1)
-v2(1)

% FR = vec2mat(FR,9);
% fr = mean(FR);
% frSTD = std(FR);
% Fr = 0.6:0.1:1.4;
% v1 = polyfit(log10(Fr(1:5)),log10(fr(1:5)),1);
% x1 = Fr(1:7);
% y1 = 10^v1(2)*x1.^v1(1);
% v2 = polyfit(log10(Fr(5:end)),log10(fr(5:end)),1);
% x2 = Fr(3:end);
% y2 = 10^v2(2)*x2.^v2(1);
% errorbar(Fr,fr,frSTD,'o')
% set(gca,'YScale','log')
% hold on
% semilogy(x1,y1,'LineWidth',1.5)
% hold on
% semilogy(x2,y2,'LineWidth',1.5)
% xlabel('Scaled IE ratio')
% ylabel('Firing Rate(Hz)')
% axes('Position',[.5 .5 .35 .35])
% box on
% errorbar(Fr,fr,frSTD,'o')
% hold on
% plot(x1,y1,'LineWidth',1.5)
% hold on
% plot(x2,y2,'LineWidth',1.5)
% -v1(1)
% -v2(1)
% % errorbar(Fr,fr,frSTD,'o-')
% xlabel('Scaled IE ratio')
% ylabel('Firing Rate(Hz)')
% end