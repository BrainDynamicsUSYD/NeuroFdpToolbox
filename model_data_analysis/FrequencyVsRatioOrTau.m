function FrequencyVsRatioOrTau
% read and converge multiple .mat to one
% specifically working in directory change_gbalance

dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
for i = 1:num_files
    R = load(files{i});
    max_lag = 200; % ms
    fre = zeros(1,16);
    for j = 1:16
        LFP_broad = R.LFP.LFP_broad(j,:);
        [ac,lags] = autocorr(LFP_broad, round(max_lag/R.dt) );
        [~,locs] = findpeaks(ac,lags*R.dt);
        if ~isempty(locs)
            fre(j) = 1e3/locs(1);
        else
            fre(j) = NaN;
        end
    end
    FRE1(i) = max(fre); % max
    [f,simple] = PowerSpectrumSpikes(R);
%     if mod(i,13) > 0 && mod(i,13) < 5
%         start = 400;
%     else
%         start = 100;
%     end
    FRE2(i) = f(simple == max(simple(30:end))); % [3,500] Hz range  % start
end
FRE1 = vec2mat(FRE1,13);
err1 = nanstd(FRE1);
FRE1 = nanmean(FRE1); % round(nanmean(FRE1))
FRE2 = vec2mat(FRE2,13);
err2 = nanstd(FRE2);
FRE2 = nanmean(FRE2);

FRE = [FRE1;FRE2];
ERR = [err1;err2];
fre = max(FRE);
err = ERR(FRE==fre);

x = 0.7:0.05:1.3;
% x = 3:0.5:9;
figure
errorbar(x,fre,err)
xlabel('F_r')
% xlabel('Decay_{GABA}(ms)')
ylabel('Frequency(Hz)')
end