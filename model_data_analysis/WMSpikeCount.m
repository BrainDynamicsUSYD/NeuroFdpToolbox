% function WMSpikeCount
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 500; % 40; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
C = zeros(1,num_files);
NumP = 3;
repeat = 100;
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
LoalNeu = cell(1,NumP);
R = load(files{1});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
SpikeC = [];
for i = 11 % 1:21
    for id_out = 1:repeat
        for no = 1:NumP
            fprintf('Processing output file No.%d out of %d...\n', i+21*(id_out-1), num_files);
            fprintf('\t File name: %s\n', files{i+21*(id_out-1)});
            R = load(files{i+21*(id_out-1)});
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},2.31e4:end)),bin,2));
            SpikeC = [SpikeC r];
        end
    end
end
%%
Signal = SpikeC(SpikeC>120); 

[N,edges] = histcounts(Signal,50);
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o')
hold on;
Y = Y(N > 0);
N = N(N > 0);
YY = Y; %(5:40);
NN = N;% (5:40);
v = polyfit(log10(YY),log10(NN),1);
x = Y; % (2:end);
y = 10^v(2)*x.^v(1);
loglog(x,y,'LineWidth',1.5)
str = ['\alpha = ',num2str(-v(1))];
text(min(x),2*max(y),str)
xlabel('Spike Count')
ylabel('Probability')
% end