function FanoFactor
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% bin = 1e3; % 100ms
%%
dir_strut2 = dir('*_config_data.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end
fprintf('\t File name: %s\n', files2{1}); % id_out-11
load(files2{1},'StiNeu')
%% fano factor of one simulation
id_out = 1;
fprintf('\t File name: %s\n', files{id_out});
R = load(files{id_out});
num = floor((4.25e4-2.25e4)/bin); % length(R.num_spikes{1}) % 4.25e4
SpikeMat = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu{1},2.25e4+1:2.25e4+bin*num)),length(StiNeu{1}),bin,[]),2))';
% num = floor(1.5e4/bin); % length(R.num_spikes{1}) % 4.25e4
% SpikeMat = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu{1},0.45e4+1:0.45e4+bin*num)),length(StiNeu{1}),bin,[]),2))';
fanoMat = var(SpikeMat)./nanmean(SpikeMat);
fanoMean = nanmean(fanoMat);
%% fano factor of multiple simulations
SpikeMat = zeros(num_files,length(StiNeu{1}));
for id_out = 1:num_files
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    SpikeMat(id_out,:) = sum(R.spike_hist{1}(StiNeu{1},3e4+1:3e4+bin),2);
    %     SpikeMat(id_out,:) = sum(R.spike_hist{1}(StiNeu{1},1e4+1:1e4+bin),2);
end
fanoMat = var(SpikeMat)./nanmean(SpikeMat);
fanoMean = nanmean(fanoMat);
%%
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
for NumP = 1:7
    switch NumP
        case 1
            Coor = [0;0];
        case 2
            Coor = [-16 15.5;-16 15.5];
        case 3
            Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
        case 4
            Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
        case 5
            Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
        case 6
            Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
                -19.6 -9.8         -9.8        9.8          9.8         19.6];
        case 7
            Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
                -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
    end
    hw = 31;
    [Lattice, ~] = lattice_nD(2, hw);
    R = load(files{1+100*(NumP-1)});
    LoalNeu = cell(1,R.ExplVar.NumP);
    for i = 1:R.ExplVar.NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    bin = 1e3:1e3:2e4;
    SpikeMat = zeros(100,length(bin));
    for i = 1:length(bin)
        for id_out = 1:100
            fprintf('Processing output file No.%d out of %d...\n', id_out+100*(NumP-1), num_files);
            fprintf('\t File name: %s\n', files{id_out+100*(NumP-1)});
            R = load(files{id_out+100*(NumP-1)});
            for no = 1:NumP
                SpikeMat(id_out,i) = SpikeMat(id_out,i) + mean(sum(R.spike_hist{1}(LoalNeu{no},3e4+1:3e4+bin(i)),2));
            end
            %     SpikeMat(id_out,:) = sum(R.spike_hist{1}(StiNeu{1},1e4+1:1e4+bin),2);
        end
    end
    fanoMat = var(SpikeMat)./nanmean(SpikeMat);
    % fanoMean = nanmean(fanoMat);
    plot(0.1*bin,fanoMat,'o-','color',Color(NumP,:))
    hold on
end
xlabel('Counting Time T (ms)')
ylabel('F(T)')
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item')
%% plotting
plot(1:7,FFinter,'bo-',1:7,FFbetween,'r>-',1:7,FFinterBefore,'g.-',1:7,FFbetweenBefore,'c.-')
legend('within one trial','across trials','within one trial before loading','across trials before loading')
xlabel('Loading Items')
ylabel('Fano Factor')
end