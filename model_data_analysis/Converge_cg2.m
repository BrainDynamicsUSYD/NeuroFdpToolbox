function Converge_cg2
% read and converge multiple .mat to one
% specifically working in directory change_gbalance
dir_strut = dir('DCAmplitudeSAmplitudePmin15msInt9000msGap1000ms-0*.mat');% -0*.mat'); % Int8000msGap1000ms-0*.mat');
% dir_strut2 = dir('AmpPatternLFPfs1000-0*.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
% files2 = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
    %     files2{id_out} = dir_strut2(id_out).name;
end
% dir_strut2 = dir('3DBurst0*minTime15SR1000P95.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
nte = zeros(1,num_files);
for i = 1:num_files
    R = load(files{i});
%     fprintf('Loading RYG.mat file %s...\n', files2{i});
%     P = load(files2{i});
%     np = length(P.ts);
%     for k = 1:np
%         send = sort(P.Send{k},'descend');
%         P.Send{k} = send(1:6);
%     end    
%     p = 1:1e3:(1e4-9e3+1);
%     for per = 1:length(p)
%         cand = find(P.ts>=p(per) & P.ts+P.Duration<=p(per)+9e3-1);
%         R.DC{per} = length(R.DC{per})/(length(cand)-1);
%     end
%     nte(i) = nanmean([R.DC{:}]);
    %     P = load(files2{i});
    %     nte(i) = length(R.dc);
    nte(i) = mean(cellfun('length',R.DC));
%         nte(i) = nanmean(cellfun(@mean,R.DC));
end

% index = vec2mat(nte,13);
% err = std(index);
% index = nanmean(index);

x = 0.7:0.05:1.3;
figure
plot(x,nte,'*-')
% errorbar(x,index,err)
title('amplitude+amplitude pattern,7E')
xlabel('IE ratio')
ylabel('Distributed Commucation Index')
end