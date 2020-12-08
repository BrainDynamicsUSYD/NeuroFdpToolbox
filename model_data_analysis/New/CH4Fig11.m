% CH4Fig11
figure_width = 14.1; % cm
figure_hight = 5.7; % cm
figure('NumberTitle','off','name', 'CH4Fig11V2', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
NumP = 3;
IE = 0.7:0.03:1.3;
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
LoalNeu = cell(1,NumP);
R = load(files{1});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
for sub = 1:2
    subplot(1,3,sub)
    MovVa = zeros(100,length(R.spike_hist{1}));
    MovVa2 = zeros(100,length(R.spike_hist{1}));
    switch sub
        case 1
            ind = 11:21:num_files;
        case 2
            ind = 1:21:num_files;
        case 3
            ind = 21:21:num_files;
    end
    for id_out = ind % 100*(NumP-1)+1:100*NumP %  % 11:21:num_files
        fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
        fprintf('\t File name: %s\n', files{id_out});
        R = load(files{id_out});
        rmat = zeros(NumP,length(R.spike_hist{1}));
        for i = 1:NumP
            rmat(i,:) = sum(movsum(full(R.spike_hist{1}(LoalNeu{i},:)),bin,2))/bin*1e4; %
            r = rmat(1,:);
        end
        MovVa(id_out,:) = r;
        MovVa2(id_out,:) = mean(rmat,1);
    end
    Fr1 = mean(MovVa(:,2.3e4:1e5),1);
    Fr1 = smooth(Fr1,1e4);
    Fr2 = mean(MovVa2(:,2.3e4:1e5),1);
    Fr2 = smooth(Fr2,1e4);
    plot((2.3e4:1e5)*1e-4,Fr1,(2.3e4:1e5)*1e-4,Fr2) %
    if sub == 1
        v2 = polyfit(log10((3e4:1e5)*1e-4),log10(Fr2(end-length(3e4:1e5)+1:end))',1);
        x2 = (2.5e4:1e5)*1e-4;
        y2 = 10^v2(2)*x2.^v2(1);
        hold on
        plot(x2,y2,'--k','LineWidth',1.5)
    end
    xlabel('Time(s)','fontSize',10)
    ylabel('Population Firng rate(Hz)','fontSize',10)
    xlim([2 9])
    switch sub
        case 1
            legend('1-item','mean','fit')
            text(-0.1,1,'A','Units', 'Normalized','FontSize',12)
        case 2
            legend('1-item','mean')
            text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
        case 3
            legend('1-item','mean')
            text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
    end
end
%%
heix = zeros(1,length(IE));
heiy = zeros(1,length(IE));
for i = 1:length(IE)
    edges = -1.5:0.05:1.5;
    [~,edges] = histcounts(Cenx{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    try
        pd = fitdist(Cenx{i}(:),'Normal'); % 'Normal'
%         y = pdf(pd,edges);
%         heix(i) = max(y/sum(y));
        heix(i) = 1/pd.sigma^2;
    catch
        heix(i) = 0;
    end
    edges = -1.5:0.05:1.5;
    [~,edges] = histcounts(Ceny{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    try
        pd = fitdist(Ceny{i}(:),'Normal');
%         y = pdf(pd,edges);
%         heiy(i) = max(y/sum(y));
        heiy(i) = 1/pd.sigma^2;
    catch
        heiy(i) = 0;
    end
end
subplot(1,3,3)
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
bin = 500; % 4ms
Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
LoalNeu = cell(1,3);
R = load(files{1});
for i = 1:R.ExplVar.NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
for r = [0.001 0.002 0.005 0.01 0.02]
    En = zeros(1,num_files);
    for id_out = 1:num_files
        fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
        fprintf('\t File name: %s\n', files{id_out});
        R = load(files{id_out});
        En(id_out) = (sum(sum(R.spike_hist{1}([LoalNeu{1}],2.26e4:end)))+r*length([LoalNeu{1}])*7.74e4/bin); % /R.ExplVar.NumP
    end
    En = vec2mat(En,length(IE));
    ee = mean(En);
    EnEf = ((heix+heiy)/2)./ee;
    plot(IE(1:end),EnEf(1:end),'o-')
    hold on
end
legend('r=0.001','r=0.002','r=0.005','r=0.01','r=0.02')
xlabel('IE ratio','fontSize',10)
ylabel('Energy Efficiency','fontSize',10)
text(-0.1,1,'C','Units', 'Normalized','fontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig11V2