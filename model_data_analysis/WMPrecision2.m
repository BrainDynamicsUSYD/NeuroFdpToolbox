% function WMPrecision2
spikec = cell(1,7);
spikecmat = cell(1,7);
spikeco = zeros(1,7);
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
post_dist = lattice_nD_find_dist(Lattice,hw,3938);
[~,IndP] = sort(post_dist);
other = IndP(1:50)';
figure
for i = 1:7
    if i > 1
        PathName = sprintf('../%dItemsSimultaneousV2',i);
        cd(PathName)
    end
    dir_strut = dir('*_RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
    dir_strut2 = dir('*_config_data.mat');
    files2 = cell(1,num_files);
    for id_out = 1:num_files
        files2{id_out} = dir_strut2(id_out).name;
    end
    load(files2{1},'StiNeu')
    for n = 1 % :length(StiNeu)
        SpikeC = zeros(1,num_files);
        SpikeCO = zeros(1,num_files);
        for j = 1:num_files
            R = load(files{j});
            SpikeC(j) = sum(sum(R.spike_hist{1}(StiNeu{n},2.26e4:7.26e4)));
            SpikeCO(j) = sum(sum(R.spike_hist{1}(other,2.26e4:7.26e4)));
        end
        spikec{i} = [spikec{i} mean(SpikeC)];
        spikecmat{i} = [spikecmat{i};SpikeC];
        spikeco(i) = mean(SpikeCO);
    end
    plot(i*ones(1,i),spikec{i},'ro')
    hold on
end
% plot(1:7,spikeco,'b>')
%%
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
spikecmat = cell(1,7);
figure
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
    LoalNeu = cell(1,NumP);
    R = load(files{100*(NumP-1)+1});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for n = 1%:NumP
        SpikeC = zeros(1,100);
        for j = 1:100
            R = load(files{100*(NumP-1)+j});
            SpikeC(j) = sum(sum(R.spike_hist{1}(LoalNeu{n},2.26e4:7.26e4)));
        end
        spikecmat{i} = [spikecmat{i};SpikeC];
    end
end
%%
figure
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
for i = 1:7 % :7
    edges = 0:50:1500;
    [N,edges] = histcounts(spikecmat{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    plot(edges,N,'o','color',Color(i,:))
    hold on
end
%
for i = 1:7 % :7
    edges = 0:50:1500;
    [N,edges] = histcounts(spikecmat{i}(:),edges,'Normalization','probability');
    edges = (edges(1:end-1)+edges(2:end))/2;
    pd = fitdist(spikecmat{i}(:),'gamma')
    %     y = gampdf(edges,pd.a,pd.b);
    %     pd = fitdist(spikecmat{i}(:),'lognormal')
    y = pdf(pd,edges);
    plot(edges,y/sum(y),'-','color',Color(i,:));
    hold on
end
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item') %
xlabel('PrecisionIndex2')
ylabel('Probability')
% title(sprintf('Load %d-Item',i))
%%
x = 1:7;
y = cellfun(@mean,spikec);
figure
v1 = polyfit(log10(x),log10(y),1);
x1 = 1:0.1:8;
y1 = 10^v1(2)*x1.^v1(1);
plot(x1,y1,'b-') %
hold on
plot(1:7,spikeco,'b>')
hold on
for i = 1:7
    plot(i*ones(1,i),spikec{i},'ro')
    hold on
end
legend('fit','Control','Data')
xlabel('Loading Items')
ylabel('Precision Index2')
% end