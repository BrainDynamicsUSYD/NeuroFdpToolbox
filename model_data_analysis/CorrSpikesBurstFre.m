function [Fre,Spikes] = CorrSpikesBurstFre
% Correlation between number of spikes in pattern with burst frequency
% Using image processing to locate bursts, which is a local method to
%  identify bursts.
tic;
dir_strut = dir('*RYG.mat');
dir_strut2 = dir('UPattern*.mat');
temp = load('0001-201804161619-42866_config_data.mat', 'LFP_centre_x');
Elecx = temp.LFP_centre_x;
temp = load('0001-201804161619-42866_config_data.mat', 'LFP_centre_y');
Elecy = temp.LFP_centre_y;
num_files = length(dir_strut);
files = cell(1,num_files);
files2 = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
    files2{id_out} = dir_strut2(id_out).name;
end
l = 63;
Fre = [];
Spikes = [];
for i = 1 %:num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    P =  load(files2{i});
    [Dfre,~,Start,End] = BurstDrifting(R);
    fprintf('Finishing burst drifting...\n')
    for j = 1:length(End) % per LFP electrode
        for k = 1:length(Start{j}) % per burst period
            Num = [];
            ind = find(P.ts >= Start{j}(k) & P.ts <= End{j}(k));
            count = length(ind);
            s = 1;
            while s < count
                d = Distance_xy(Elecx,Elecy,P.x(ind(s)),P.y(ind(s)),l);
                if min(d) == d(j)
                    Num = [Num P.num(ind(s))];
                end
                s = s + 1;
            end
            if ~isempty(Num)
                Fre = [Fre Dfre{j}(k)];
                Spikes = [Spikes mean(Num)];
            end
            toc;
        end
        toc;
    end
end
toc;
end