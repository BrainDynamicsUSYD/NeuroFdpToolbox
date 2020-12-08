function SynapticIXCorr
dir_strut = dir('*_0_neurosamp.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
Lag = 500; % 50 ms
for i = 4 % 1:num_files
    load(files{i},'I_AMPA','I_GABA');
    X = zeros(size(I_AMPA,1),2*Lag+1);
    L = zeros(size(I_AMPA,1),2*Lag+1);
    for j = 1:size(I_AMPA,1)
        [xcor,lag] = xcorr(I_AMPA(j,:),-I_GABA(j,:),Lag);
        X(j,:) = xcor;
        L(j,:) = lag;
    end
    plot(0.1*mean(L),mean(X))
    xlabel('Lag Time(ms)')
    ylabel('XCorrelation')
    %         next = input('\t Next figure?');
    %         delete(gcf);
end
end