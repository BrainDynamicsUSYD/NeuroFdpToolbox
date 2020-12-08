function WMEnergyEfficiency
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
% bin = 40; % 4ms
%% mean(binary)/energy
% EnEf = zeros(1,num_files); % num_files
% for id_out = 1:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
%     load(files2{id_out},'StiNeu')
%     c = 0;
%     for i = 1:length(StiNeu)
%         r = sum(movsum(full(R.spike_hist{1}(StiNeu{i},2.25e4+1:end)),bin,2));
%         if max(r) > 25
%             c = c + 1;
%         end
%     end
%     EnEf(id_out) = c/sum(R.num_spikes{1}); % (2.25e4+1:end)
% %     EnEf(id_out) = c/sum(sum(R.spike_hist{1}(StiNeu{i},:)));
% % a = [StiNeu{:}];
% % EnEf(id_out) = c/sum(sum(R.spike_hist{1}(a(:),:)));
% end
% EnEf = vec2mat(EnEf,9);
% ee = mean(EnEf);
% ee_std = std(EnEf);
% x = 0.6:0.1:1.4;
% %
% errorbar(x,ee,ee_std,'o-')
% % plot(x,EnEf,'o-')
% xlim([0.5 1.5])
% xlabel('IE balance')
% ylabel('Efficiency')
% end

%% mean(binary)*height_precision/energy
% Cen = cell(1,num_files); % num_files
% c = zeros(1,num_files);
% En = zeros(1,num_files);
% hw = 31;
% [Lattice,~] = lattice_nD(2, hw);
% no = 1;
% edges = linspace(-1,1,31); % 1:[-1,1]; 2:[-24,-22]; 3&4&5:[-17,-15];
% for id_out = 1:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
%     load(files2{id_out},'StiNeu')
%     RecallC = [];
%     r = sum(movsum(full(R.spike_hist{1}(StiNeu{no},:)),bin,2));
%     candx = Lattice(StiNeu{no},1);
%     candy = Lattice(StiNeu{no},2);
%     for i = 2.76e4:length(r)-bin/2 % 4.46e4 % length(r) %
%         if r(i) > 25
%             % neurons within WM area
%             spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
%             spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
%             spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
%             spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
%
%             [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
%             RecallC = [RecallC;x y length(spike_x_pos_o) i];
%         end
%     end
%     % Recall Centers
%     xr = [];
%     yr = [];
%     xtemp = [];
%     ytemp = [];
%     num = [];
%     for i = 1:size(RecallC,1)
%         xtemp = [xtemp RecallC(i,1)];
%         ytemp = [ytemp RecallC(i,2)];
%         num = [num RecallC(i,3)];
%         if (i==size(RecallC,1)) || RecallC(i+1,4)-RecallC(i,4) > 1
%             ind = find(num == max(num));
%             xr = [xr mean(xtemp(ind))];
%             yr = [yr mean(ytemp(ind))];
%             xtemp = [];
%             ytemp = [];
%             num = [];
%         end
%     end
%     Cen{id_out} = [xr;yr];
%     En(id_out) = sum(R.num_spikes{1});
%     c(id_out) = ~isempty(Cen{id_out});
% end
% Pr = zeros(1,9);
% for i = 1:9
%     r = [Cen{10*i-9:10*i}];
%     [xN,edges] = histcounts(r(1,:),edges,'Normalization','probability');
%     [yN,edges] = histcounts(r(2,:),edges,'Normalization','probability');
%     p = (edges(1:end-1)+edges(2:end))/2;
%     pdx = fitdist(r(1,:)','Normal');
%     yx = pdf(pdx,p);
%     pdy = fitdist(r(2,:)','Normal');
%     yy = pdf(pdy,p);
%     Pr(i) = (max(yx/sum(yx))+max(yy/sum(yy)))/2;
% end
% en = mean(vec2mat(En,9));
% ci = sum(vec2mat(c,9))/10;
% x = 0.6:0.1:1.4;
% plot(x,Pr.*ci./en,'o-')
% xlim([0.5 1.5])
% xlabel('IE balance')
% ylabel('Efficiency')

%%
% Pr = zeros(1,9);
% for i = 1:9
%     r = [Cen{10*i-9:10*i}];
%     [xN,edges] = histcounts(r(1,:),edges,'Normalization','probability');
%     [yN,edges] = histcounts(r(2,:),edges,'Normalization','probability');
%     p = (edges(1:end-1)+edges(2:end))/2;
%     pdx = fitdist(r(1,:)','Normal');
%     yx = pdf(pdx,p);
%     pdy = fitdist(r(2,:)','Normal');
%     yy = pdf(pdy,p);
%     Pr(i) = (max(yx/sum(yx))+max(yy/sum(yy)))/2;
% end
% en = mean(vec2mat(En,9));
% ci = sum(vec2mat(c,9))/10;
% x = 0.6:0.1:1.4;
% plot(x,Pr.*ci./en,'o-')
% xlim([0.5 1.5])
% xlabel('IE balance')
% ylabel('Efficiency')
%% fr_precision/energy
% figure
% for r = [0.001 0.002 0.005 0.01 0.02]
% for i = 3%1:5
%     EnEf = zeros(1,num_files);
%     for id_out = 1:num_files
%         fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%         fprintf('\t File name: %s\n', files{id_out});
%         R = load(files{id_out});
%         fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
%         load(files2{id_out},'StiNeu')
%         EnEf(id_out) = sum(sum(R.spike_hist{1}(StiNeu{1},2.26e4:2.26e4+1e4*i)))/i/(sum(R.num_spikes{1})+3969*r*length(R.num_spikes{1})/40);
%     end
%     EnEf = vec2mat(EnEf,9);
%     ee = mean(EnEf);
%     ee_std = std(EnEf);
%     x = 0.6:0.1:1.4;
%     errorbar(x,ee,ee_std,'o-')
%     hold on
% end
% end
% legend('r=0.001','r=0.002','r=0.005','r=0.01','r=0.02')
% % legend('1 s','2 s','3 s','4 s','5 s')
% xlim([0.5 1.5])
% xlabel('IE balance')
% ylabel('Efficiency')
%%
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
IE = 0.7:0.03:1.3;
% Cenx = cell(1,length(IE));
% Ceny = cell(1,length(IE));
NumP = 3;
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
R = load(files{1});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
%%
for j = 1:length(IE)
    for no = 1:NumP
        for id_out = 1:100
            fprintf('Processing output file No.%d out of %d...\n', length(IE)*(id_out-1)+j, num_files);
            fprintf('\t File name: %s\n', files{length(IE)*(id_out-1)+j});
            R = load(files{length(IE)*(id_out-1)+j});
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
            candx = Lattice(LoalNeu{no},1);
            candy = Lattice(LoalNeu{no},2);
            for i = 2.26e4:length(r)-bin/2 % ind % 4.46e4 % length(r) %
                if r(i) >= 120 % 120 % floor(0.5*length(candx)) %
                    % neurons within WM area
                    spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
                    Cenx{j} = [Cenx{j} x-Coor(1,no)];
                    Ceny{j} = [Ceny{j} y-Coor(2,no)];
                end
            end
        end
    end
end
save('CenxyMat.mat','Cenx','Ceny')
%% precision distribution
% heix = zeros(1,length(IE));
figure
IE = 1:7;
n = 50;
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
for i = 1:length(IE)
    heiy = zeros(1,n);
    %     edges = -1.5:0.05:1.5;
    %     [~,edges] = histcounts(Cenx{i}(:),edges,'Normalization','probability'); % ,edges
    %     edges = (edges(1:end-1)+edges(2:end))/2;
    %     try
    %         pd = fitdist(Cenx{i}(:),'Normal'); % 'Normal'
    % %         y = pdf(pd,edges);
    % %         heix(i) = max(y/sum(y));
    %         heix(i) = 1/pd.sigma^2;
    %     catch
    %         heix(i) = 0;
    %     end
    int = floor(length(Ceny{i})/n);
    for j = 1:n
        %     try
        pd = fitdist(Ceny{i}((j-1)*int+1:j*int)','Normal');
        heiy(j) = 1/pd.sigma^2;
        %     catch
        %         heiy(j) = 0;
        %     end
    end
    edges = 0:4:200;
    [N,edges] = histcounts(heiy,edges,'Normalization','probability');
    edges = (edges(1:end-1)+edges(2:end))/2;
    plot(edges,N,'o','color',Color(i,:))
    hold on
end
for i = 1:length(IE)
    heiy = zeros(1,n);
    int = floor(length(Ceny{i})/n);
    for j = 1:n
        %     try
        pd = fitdist(Ceny{i}((j-1)*int+1:j*int)','Normal');
        heiy(j) = 1/pd.sigma^2;
        %     catch
        %         heiy(j) = 0;
        %     end
    end
    edges = 0:4:200;
    [N,edges] = histcounts(heiy,edges,'Normalization','probability');
    edges = (edges(1:end-1)+edges(2:end))/2;
    pd = fitdist(heiy','gamma')
    y = pdf(pd,edges);
    plot(edges,y/sum(y),'-','color',Color(i,:));
    hold on
end
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item') %
xlabel('Precision')
ylabel('Probability')
%% precision distribution2
figure
IE = 1:7;
n = 50;
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
for i = 1:length(IE)
    heiy = zeros(1,n);
    for j = 1:n
        try
        pd = fitdist([Ceny{i,round(1000/n*(j-1)+1:1000/n*j)}]','Normal');
        heiy(j) = 1/pd.sigma^2;
        catch
            heiy(j) = 0;
        end
    end
    edges = 0:4:200;
    [N,edges] = histcounts(heiy,edges,'Normalization','probability');
    edges = (edges(1:end-1)+edges(2:end))/2;
    plot(edges,N,'o','color',Color(i,:))
    hold on
end
for i = 1:length(IE)
    heiy = zeros(1,n);
    for j = 1:n
        try
        pd = fitdist([Ceny{i,round(1000/n*(j-1)+1:1000/n*j)}]','Normal');
        heiy(j) = 1/pd.sigma^2;
        catch
            heiy(j) = 0;
        end
    end
    edges = 0:4:200;
    [N,edges] = histcounts(heiy,edges,'Normalization','probability');
    edges = (edges(1:end-1)+edges(2:end))/2;
    pd = fitdist(heiy','gamma')
    y = pdf(pd,edges);
    plot(edges,y/sum(y),'-','color',Color(i,:));
    hold on
end
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item') %
xlabel('Precision')
ylabel('Probability')
%% heix/y: height or 1/STD^2
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
%% normalized precision fit
NP = heiy; % heiy; % (heix + heiy)/2;
NP = NP/NP(1);
v2 = polyfit(log10(IE),log10(NP),1);
x2 = IE;
y2 = 10^v2(2)*x2.^v2(1);
plot(IE,NP,'o',x2,y2,'LineWidth',1.5)
xlabel('WM items')
ylabel('Relative Precision')
legend('precision','power-law fit')
%%
figure
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
bin = 500; % 4ms
for r = [0.001 0.002 0.005 0.01 0.02]
    En = zeros(1,num_files);
    for id_out = 1:num_files
        fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
        fprintf('\t File name: %s\n', files{id_out});
        R = load(files{id_out});
        if mod(id_out,1000) == 1
            switch 7 % R.ExplVar.NumP
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
            LoalNeu = cell(1,R.ExplVar.NumP);
            R = load(files{id_out});
            for i = 1:R.ExplVar.NumP
                dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
                LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
            end
        end
        En(id_out) = (sum(sum(R.spike_hist{1}([LoalNeu{:}],2.26e4:end)))+r*length([LoalNeu{:}])*7.74e4/bin)/R.ExplVar.NumP;
    end
%     En = vec2mat(En,length(IE));
%     ee = mean(En);
    En = vec2mat(En,1000);
    ee = mean(En,2)';
    EnEf = ((heix+heiy)/2)./ee;
    plot(IE(1:end),EnEf(1:end),'o-')
    hold on
end
legend('r=0.001','r=0.002','r=0.005','r=0.01','r=0.02')
xlabel('WM items')
ylabel('Efficiency')
% saveas(gcf,'WMEnergyEfficiencyAll.eps')
end