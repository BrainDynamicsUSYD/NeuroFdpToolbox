% function WMVariability
% Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
% bin = 40; % 4ms
% % shift = [0 23 16*ones(1,5)];
% % shift2 = [0 23 16 16 16 8 8];
% for ifolder = 1:7
%     if ifolder > 1
%         PathName = sprintf('../%dItemsSimultaneousV2',ifolder);
%         cd(PathName)
%     end
%     dir_strut = dir('*_RYG.mat');
%     num_files = length(dir_strut);
%     files = cell(1,num_files);
%     for id_out = 1:num_files
%         files{id_out} = dir_strut(id_out).name;
%     end
%     dir_strut2 = dir('*_config_data.mat');
%     files2 = cell(1,num_files);
%     for id_out = 1:num_files
%         files2{id_out} = dir_strut2(id_out).name;
%     end
%     load(files2{1},'StiNeu')
%     Cen = zeros(2,num_files); % num_files
%     hw = 31;
%     [Lattice,~] = lattice_nD(2, hw);
%     for id_out = 1:num_files
%         fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%         fprintf('\t File name: %s\n', files{id_out});
%         R = load(files{id_out});
%         fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
%         load(files2{id_out},'StiNeu','IndC')
%         spike_x = [];
%         spike_y = [];
%         for no = 1:length(StiNeu)
%             %             r = sum(movsum(full(R.spike_hist{1}(StiNeu{no},:)),bin,2));
%             candx = Lattice(StiNeu{no},1);
%             candy = Lattice(StiNeu{no},2);
%             %             ind = find(r > 25);
%             for i = 2.26e4:5.26e4 % length(R.num_spikes{1}) % ind % 4.46e4 % length(r) %
%                 if 1% i >= 2.76e4
%                     % neurons within WM area
%                     spike_x_pos_o = candx.*R.spike_hist{1}(StiNeu{no},i);
%                     spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(StiNeu{no},i));
%                     spike_y_pos_o = candy.*R.spike_hist{1}(StiNeu{no},i);
%                     %                     spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(StiNeu{1},i);
%                     spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(StiNeu{no},i));
%                     spike_x = [spike_x;spike_x_pos_o-Lattice(IndC(no),1)];
%                     spike_y = [spike_y;spike_y_pos_o-Lattice(IndC(no),2)];
%                 end
%             end
%         end
%         try
%             [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x,spike_y,2*hw+1,'quick');
%         catch
%             x = NaN;
%             y = NaN;
%         end
%         Cen(:,id_out) = [x;y];
%     end
%
%     %     Cm = min(Cen(:));
%     %     Cen = Cen + abs(Cm) + 1;
%     %     Cen = Cen/max(Cen(:));
%
%     %     for c = 2 %:2
%     %         [N,edges] = histcounts(Cen(c,:),20,'Normalization','probability');
%     %         edges = (edges(1:end-1)+edges(2:end))/2;
%     %         figure(1)
%     %         subplot(2,4,ifolder)
%     %         plot(edges,N,'o-','color',Color(ifolder,:));
%     %         TitleName = sprintf('%d-item',ifolder);
%     %         title(TitleName)
%     %         figure(2)
%     %         hold on
%     %         plot(edges+shift2(ifolder),N,'o-','color',Color(ifolder,:));
%     % %         pd = fitdist(Cen(c,:)','stable')
%     % %         y = pdf(pd,edges);
%     % %         plot(edges-(abs(Cm) + 1),y/sum(y),'-','color',Color(ifolder,:));
%     % %         hold on
%     %     end
%     d = Distance_xy(Cen(1,:),Cen(2,:),0,0,63);
%     [N,edges] = histcounts(d,20,'Normalization','probability');
%     edges = (edges(1:end-1)+edges(2:end))/2;
%     plot(edges,N,'o','color',Color(ifolder,:)); % ,'MarkerSize',12
%     hold on
%     pd = fitdist(d','gamma')
%     p = pdf(pd,edges);
%     plot(edges,p/sum(p),'-','color',Color(ifolder,:));
%     hold on
% end
% xlabel('Error')
% ylabel('Probability')
% legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item')
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
Cenx = cell(1,7);
Ceny = cell(1,7);
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
    R = load(files{1+100*(NumP-1)});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for no = 1%:NumP
        for id_out = 1:100
            fprintf('Processing output file No.%d out of %d...\n', id_out+100*(NumP-1), num_files);
            fprintf('\t File name: %s\n', files{id_out+100*(NumP-1)});
            R = load(files{id_out+100*(NumP-1)});
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
%             if flag == 0
%                 spike_x = [];
%                 spike_y = [];
%             end
            candx = Lattice(LoalNeu{no},1);
            candy = Lattice(LoalNeu{no},2);
            for i = 2.26e4:length(r)-bin/2 % ind % 4.46e4 % length(r) %
                if r(i) >= 120 % floor(0.5*length(candx)) % 
                    % neurons within WM area
                    spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
                    Cenx{NumP} = [Cenx{NumP} x-Coor(1,no)];
                    Ceny{NumP} = [Ceny{NumP} y-Coor(2,no)];
%                     spike_x = [spike_x x-Coor(1,no)];
%                     spike_y = [spike_y y-Coor(2,no)];
                    %                     spike_x_pos_o = candx.*R.spike_hist{1}(LoalNeu{no},i);
                    %                     spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i));
                    %                     spike_y_pos_o = candy.*R.spike_hist{1}(LoalNeu{no},i);
                    %                     spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i));
                    %                     spike_x = [spike_x;spike_x_pos_o-Coor(1,no)];
                    %                     spike_y = [spike_y;spike_y_pos_o-Coor(2,no)];
                end
            end
%             if mod(id_out,2) == 0
%                 try
%                     pd = fitdist(spike_x(:),'Normal');
%                     Cenx{NumP} = [Cenx{NumP} pd.sigma];
%                 catch
%                 end
%                 try
%                     pd = fitdist(spike_y(:),'Normal');
%                     Ceny{NumP} = [Ceny{NumP} pd.sigma];
%                 catch
%                 end
%                 flag = 0;
%             else
%                 flag = 1;
%             end
            
            %             try
            %                 [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x,spike_y,2*hw+1,'quick');
            %             catch
            %                 x = NaN;
            %                 y = NaN;
            %             end
            %             Cen(:,id_out) = [x;y];
        end
    end
    %     d = Distance_xy(Cen(1,:),Cen(2,:),0,0,63);
    %     [N,edges] = histcounts(d,20,'Normalization','probability');
    %     edges = (edges(1:end-1)+edges(2:end))/2;
    %     plot(edges,N,'o','color',Color(ifolder,:)); % ,'MarkerSize',12
    %     hold on
    %     pd = fitdist(d','gamma')
    %     p = pdf(pd,edges);
    %     plot(edges,p/sum(p),'-','color',Color(ifolder,:));
    %     hold on
end
%%
ff = zeros(1,7);
for i = 1:7 % :7
    ff(i) = var(1./Cenx{i}(:).^2)./mean(1./Cenx{i}(:).^2);
end
plot(1:7,ff,'o-')
xlabel('WM items')
ylabel('Fano factor')
%%
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
for i = 1:7 % :7
    edges = -1.5:0.05:1.5;
    [N,edges] = histcounts(Cenx{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    plot(edges,N,'o','color',Color(i,:))
    hold on
end
%
heix = zeros(1,7);
STDx = zeros(1,7);
for i = 1:7 % :7
    edges = -1.5:0.05:1.5;
    [N,edges] = histcounts(Cenx{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    pd = fitdist(Cenx{i}(:),'Normal') % 'Normal'
    STDx(i) = pd.sigma;
    y = pdf(pd,edges);
    heix(i) = max(y/sum(y));
    plot(edges,y/sum(y),'-','color',Color(i,:));
    hold on
end
xlabel('Errorx')
ylabel('Probability')
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item')
%
figure
plot(1:7,STDx,'o-b',1:7,heix,'*-g')
xlabel('WM items')
ylabel('ErrorxIndex')
legend('STD','Height')
figure
for i = 1:7 % :7
    edges = -1.5:0.05:1.5;
    [N,edges] = histcounts(Ceny{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    plot(edges,N,'o','color',Color(i,:))
    hold on
end
%
heiy = zeros(1,7);
STDy = zeros(1,7);
for i = 1:7 % :7
    edges = -1.5:0.05:1.5;
    [N,edges] = histcounts(Ceny{i}(:),edges,'Normalization','probability'); % ,edges
    edges = (edges(1:end-1)+edges(2:end))/2;
    pd = fitdist(Ceny{i}(:),'Normal')
    STDy(i) = pd.sigma;
    y = pdf(pd,edges);
    heiy(i) = max(y/sum(y));
    plot(edges,y/sum(y),'-','color',Color(i,:));
    hold on
end
xlabel('Errory')
ylabel('Probability')
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item')
figure
plot(1:7,STDy,'o-b',1:7,heiy,'*-g')
xlabel('WM items')
ylabel('ErroryIndex')
legend('STD','Height')
%
figure
plot(1:7,(STDx+STDy)/2,'o-b',1:7,(heix+heiy)/2,'*-g')
xlabel('WM items')
ylabel('ErrorIndex')
legend('STD','Height')
%%
Color = [1 1 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1;0 0 0];
STDymat = zeros(7,10);
for i = 1:7 % :7
    group = 1:floor(length(Ceny{i}(:))/10):length(Ceny{i}(:));
    for j = 1:length(group)-1
    pd = fitdist(Ceny{i}(group(j):group(j+1))','Normal'); % 'Normal'
    STDymat(i,j) = pd.sigma;
    end
end
figure
s = std(STDymat,0,2);
errorbar(1:7,STDy,s,'o-')
xlabel('WM items')
ylabel('ErroryIndex')
%% movsum variability
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
for NumP = 1:7
% IE = 0.7:0.03:1.3;
% switch NumP
%     case 1
%         Coor = [0;0];
%     case 2
%         Coor = [-16 15.5;-16 15.5];
%     case 3
%         Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
%     case 4
%         Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
%     case 5
%         Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
%     case 6
%         Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%             -19.6 -9.8         -9.8        9.8          9.8         19.6];
%     case 7
%         Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%             -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
% end
% LoalNeu = cell(1,NumP);
R = load(files{1});
% for i = 1:NumP
%     dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
%     LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
% end
MovVa = zeros(100,length(R.spike_hist{1}));
% MovVa2 = zeros(100,length(R.spike_hist{1}));
% ind = 1:21:num_files;
for id_out = 1:100 % 100*(NumP-1)+1:100*NumP %  % 11:21:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{100*(NumP-1)+id_out});
    R = load(files{100*(NumP-1)+id_out});
    r = sum(movsum(full(R.spike_hist{1}),bin,2))/bin*1e4/3969;
%     rmat = zeros(NumP,length(R.spike_hist{1}));
%     for i = 1:NumP
%         rmat(i,:) = sum(movsum(full(R.spike_hist{1}(LoalNeu{i},:)),bin,2))/bin*1e4; % 
%         r = rmat(1,:);
%     end
    MovVa(id_out,:) = r;
%     MovVa2(id_out,:) = mean(rmat,1);
end
%
Fr1 = mean(MovVa(:,2.5e4:1e5),1);
Fr1 = smooth(Fr1,1e4);
% Fr2 = mean(MovVa2(:,2.5e4:1e5),1);
% Fr2 = smooth(Fr2,1e4);
% v2 = polyfit(log10((3e4:1e5)*1e-4),log10(Fr2(end-length(3e4:1e5)+1:end))',1);
% x2 = (2.5e4:1e5)*1e-4;
% y2 = 10^v2(2)*x2.^v2(1);
plot((2.5e4:1e5)*1e-4,Fr1,'color',Color(NumP,:)) % ,(2.5e4:1e5)*1e-4,Fr2)
hold on
end
% plot(x2,y2,'--k','LineWidth',1.5)
% title(sprintf('%s=%0.2f','\alpha',v2(1))) 
xlabel('Time(s)')
ylabel('Population Firng rate(Hz)')
legend('1-item','2-item','3-item','4-item','5-item','6-item','7-item')
xlim([2 9])
% legend('1-item','mean') % ,'fit')
% end