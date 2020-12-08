% function Converge_decay
% read and converge multiple .mat to one
% specifically working in directory change_gbalance
dir_strut = dir('DIV2*.mat');
dir_strut2 = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
files2 = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
    files2{id_out} = dir_strut2(id_out).name;
end
index1 = [];
index2 = [];
for i = 1:num_files
    R = load(files{i});
    R2 = load(files2{i});
    index1 = [index1 R.DC(1)];
    index2 = [index2 R.DC(2)];
    max_lag = 200; % ms
    fre = zeros(1,16);
    for j = 1:16
        LFP_broad = R2.LFP.LFP_broad(j,:);
        [ac,lags] = autocorr(LFP_broad, round(max_lag/R2.dt) );
        [~,locs] = findpeaks(ac,lags*R2.dt);
        if ~isempty(locs)
            fre(j) = 1e3/locs(1);
        else
            fre(j) = NaN;
        end
    end
    FRE(i) = max(fre);
end
FRE = round(FRE);
% labels = cellstr(num2str(FRE'));
ind = [1:11 13];
index1(ind) = 9/19*index1(ind);
index2(ind) = 9/19*index2(ind);
x = FRE;
% disp(x)
[x,I] = sort(x);
index1 = index1(I);
index2 = index2(I);
labels = cellstr(num2str(x'));
figure
plot(x,index1,x,index2)
text(x,index1,labels)
text(x,index2,labels)
legend('index 1','index 2')
xlabel('Frequency')
ylabel('normalized distributed commucation index')
% end