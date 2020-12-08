function WMPrecision
% using spiking number account for precision
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
dir_strut2 = dir('*_config_data.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
load(dir_strut2.name,'StiNeu','IndC')
% NumP = length(StiNeu);
% r = cell(1,NumP);
% Color = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0];

%% small bin
% bin = 400; % 4ms
% for id_out = 1:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out-11}); % id_out-11
%     load(files2{id_out-11},'StiNeu')
%     NumP = length(StiNeu);
%     %     r = cell(1,NumP);
%     for i = 1:NumP
%         %         r{i} = sum(movsum(full(R.spike_hist{1}(StiNeu{i},:)),bin,2));
%         %         subplot(NumP,1,i)
%         %         plot(0.1*(1:R.step_tot),r{i})
%         %         ylim([0 max(r{i})])
%         %         xlabel('Time(ms)')
%         %         ylabel('Spike Counts')
%
%         % average over multiple trials
%         r{i} = [r{i};sum(movsum(full(R.spike_hist{1}(StiNeu{i},:)),bin,2))];
%     end
%     %     next = input('\t Next figure?');
%     %     close all
% end
% for i = 1:NumP
%     subplot(NumP,1,i)
%     plot(0.1*(1:length(r{i})),mean(r{i}))
%     ylim([0 length(StiNeu{1})])
%     xlabel('Time(ms)')
%     ylabel('Spike Counts')
% end

%% sum of spikes--different items in single trial
% Color = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0];
% for id_out = 15:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out-7}); % id_out-11/7
%     load(files2{id_out-7},'StiNeu')
%     NumP = length(StiNeu);
%     r = cell(1,NumP);
%     for i = 1:NumP
%         r{i} = sum(full(R.spike_hist{1}(StiNeu{i},2.35e4+1:end)),2);
%         [N,edges] = histcounts(r{i});
%         edge = round((edges(1:end-1)+edges(2:end))/2);
%         plot(edge,N,'o-','color',Color(i,:))
%         hold on
%         xlabel('Spike Counts')
%         ylabel('Neuron Counts')
%     end
%     next = input('\t Next figure?');
%     close all
% end

%% sum of spikes--multiple items over single trial
% r = [];
% r2 = [];
% for id_out = 1:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     %     fprintf('\t File name: %s\n', files2{id_out});
%     %     load(files2{id_out},'StiNeu','IndC')
%     %     NumP = length(StiNeu);
%     for i = 1:NumP
%         r = [r sum(full(R.spike_hist{1}(StiNeu{i},2.25e4+1:end-5e4)),2)];
%         r2 = [r2 sum(full(R.spike_hist{1}(StiNeu{i},2.25e4+1:end)),2)];
%     end
% end
% SC = mean(mean(r)); % spike count
% SC2 = mean(mean(r2));
% for i = 1:NumP
%     [N,edges] = histcounts(mean(r{i},2));
%     edge = round((edges(1:end-1)+edges(2:end))/2);
%     plot(edge,N,'o-','color',Color(i,:))
%     hold on
%     xlabel('Spike Counts')
%     ylabel('Neuron Counts')
% end
% Neu = cat(1,StiNeu{:});
% hw = 31;
% [Lattice, ~] = lattice_nD(2, hw);
% plot3(Lattice(Neu,1),Lattice(Neu,2),r,'.')
% zlabel('Spike Count')
% grid on
%%
% figure
% dist = [];
% for i = 1:NumP
%     dist = [dist; Distance_xy(Lattice(StiNeu{i},1),Lattice(StiNeu{i},2),Lattice(IndC(i),1),Lattice(IndC(i),2),2*hw+1)];
% end
% plot(repmat(dist,num_files,1),r,'.')
% xlabel('Distance from the center')
% ylabel('Spike Count')
%% sum of spikes tuning curve on x axis
% NumP = 4;
% r = cell(1,NumP);
% r2 = cell(1,NumP);
% for id_out = 1:num_files
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out});
%     load(files2{id_out},'StiNeu','IndC')
%     for i = 1:NumP
%         r{i} = [r{i} sum(full(R.spike_hist{1}(IndC(i)-10:IndC(i)+10,2.25e4+1:end-5e4)),2)];
%         r2{i} = [r2{i} sum(full(R.spike_hist{1}(IndC(i)-10:IndC(i)+10,2.25e4+1:end)),2)];
%         %         if id_out == 100
%         %             subplot(2,2,i)
%         %             plot(-10:10,mean(r{i},2),'o-',-10:10,mean(r2{i},2),'*-')
%         %             legend('2.75s-period','7.75s-period')
%         %             ylim([0 35])
%         %             xlabel('Distance From the Center(a.u.)')
%         %             ylabel('Spike Count')
%         %         end
%     end
% end
% figure
% plot(-10:10,mean(cat(2,r{:}),2),'o-',-10:10,mean(cat(2,r2{:}),2),'*-')
% legend('2.75s-period','7.75s-period')
% x = -10:10; % -10:10 % -4:4
% y1 = mean(cat(2,r{:}),2)';
% % y1 = y1(7:15);
% y2 = mean(cat(2,r2{:}),2)';
% % y2 = y2(7:15);
% f1 = fit(x',y1','gauss1');
% f2 = fit(x',y2','gauss1');
% hold on;plot(f1,x,y1)
% hold on;plot(f2,x,y2)
% ylim([0 30])
% xlabel('Distance From the Center(a.u.)')
% ylabel('Spike Count')
%% center of spikes
bin = 40; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
xr = [];
yr = [];
no = 1;
% subplot(2,2,4)
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
%     fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
%     load(files2{id_out},'StiNeu')
    LoadC = [];
    RecallC = [];
    OtherC = [];
    r = sum(movsum(full(R.spike_hist{1}(StiNeu{no},:)),bin,2));
    candx = Lattice(StiNeu{no},1);
    candy = Lattice(StiNeu{no},2);
    for i = 2.5e4:2.75e4
        if r(i) > 30 % 30 % 60pr
            spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(StiNeu{no},i-bin/2+1:i+bin/2);
            spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(StiNeu{no},i-bin/2+1:i+bin/2));
            spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(StiNeu{no},i-bin/2+1:i+bin/2);
            spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(StiNeu{no},i-bin/2+1:i+bin/2));
            [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
            LoadC = [LoadC;x y length(spike_x_pos_o)];
        end
    end
    for i = 2.76e4:length(r)-bin/2 % 4.46e4 % length(r) % 
        if r(i) > 45 % 25
            % neurons within WM area
%             spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
%             spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
%             spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2);
%             spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(StiNeu{1},i-bin/2+1:i+bin/2));
            % all neurons
            spike_x_pos_o = repmat(Lattice(:,1),1,bin).*R.spike_hist{1}(:,i-bin/2+1:i+bin/2);
            spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(:,i-bin/2+1:i+bin/2));
            spike_y_pos_o = repmat(Lattice(:,2),1,bin).*R.spike_hist{1}(:,i-bin/2+1:i+bin/2);
            spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(:,i-bin/2+1:i+bin/2));
            
            [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
            RecallC = [RecallC;x y length(spike_x_pos_o) i];
        end
    end
    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(candx,candy,2*hw+1,'quick');
    %     plot(RecallC(:,1),RecallC(:,2),'bo',mean(LoadC(:,1)),mean(LoadC(:,2)),'r>',x,y,'r*')
    %     plot(RecallC(:,1),RecallC(:,2),'bo',xl,yl,'r>',x,y,'r*')
    % Load Center
    ind = find(LoadC(:,3) == max(LoadC(:,3)));
    xl = mean(LoadC(ind,1));
    yl = mean(LoadC(ind,2));
    % Recall Centers
    %     xr = [];
    %     yr = [];
    xtemp = [];
    ytemp = [];
    num = [];
    for i = 1:size(RecallC,1)
        xtemp = [xtemp RecallC(i,1)];
        ytemp = [ytemp RecallC(i,2)];
        num = [num RecallC(i,3)];
        if (i==size(RecallC,1)) || RecallC(i+1,4)-RecallC(i,4) > 1
            ind = find(num == max(num));
            xr = [xr mean(xtemp(ind))];
            yr = [yr mean(ytemp(ind))];
            xtemp = [];
            ytemp = [];
            num = [];
        end
    end
    %     if ~isempty(RecallC)
    %         other = setdiff(2.26e4:4.46e4,RecallC(:,4));
    %     else
    %         other = 2.26e4:4.46e4;
    %     end
    %     while length(other) > bin
    %         if other(bin)-other(1) == (bin-1)
    %             spike_x_pos_o = repmat(Lattice(:,1),1,bin).*R.spike_hist{1}(:,other(1):other(bin));
    %             spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(:,other(1):other(bin)));
    %             spike_y_pos_o = repmat(Lattice(:,2),1,bin).*R.spike_hist{1}(:,other(1):other(bin));
    %             spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(:,other(1):other(bin)));
    %             if ~isempty(spike_x_pos_o)
    %                 [xo,yo,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
    %                 OtherC = [OtherC;xo yo];
    %             end
    %             other = other(bin+1:end);
    %         else
    %             other = other(2:end);
    %         end
    %     end
    %     plot(OtherC(:,1),OtherC(:,2),'g.',xr,yr,'bo',xl,yl,'r>',x,y,'r*')
%     plot(xr,yr,'bo',xl,yl,'r>',x,y,'r*')
%     xlim([-31,31]) % ([-17 -15])
%     ylim([-31,31]) % ([-17 -15])
%     hold on;
%     next = input('\t Next figure?');
%     close all
end
plot(xr,yr,'bo',xl,yl,'r>',x,y,'r*')
xlim([-31,31]) % ([-31,31]) % ([-17 -15])
ylim([-31,31]) % ([-31,31]) % ([-17 -15])
% hold on;
% axes('Position',[.5 .5 .35 .35])
% box on
% plot(xr,yr,'bo',xl,yl,'r>',x,y,'r*')
% xlim([-17 -15]) % ([-1,1]) % ([-17 -15])
% ylim([-17 -15]) % ([-1,1]) % ([-17 -15])

%%
edges = linspace(-17,-15,31); % 1:[-1,1]; 2:[-24,-22]; 3&4&5:[-17,-15];
[xN,edges] = histcounts(xr,edges,'Normalization','probability');
[yN,edges] = histcounts(yr,edges,'Normalization','probability');
% 6&7 items
% [xN,edges] = histcounts(xr,linspace(-17,-15,31),'Normalization','probability');

p = (edges(1:end-1)+edges(2:end))/2;
subplot(1,2,1)
pd = fitdist(xr','Normal');
y = pdf(pd,p);
plot(p,xN,'b.',p,y/sum(y),'r--')
xlabel('Position on x axis')
ylabel('Probability')
% str = [{['Width=',num2str(pd.sigma)]},{['Height=',num2str(max(xN))]},{['FitHeight=',num2str(max(y/sum(y)))]}];
str = sprintf('Width=%.2f\nHeight=%.2f\nFitHeight=%.2f',pd.sigma,max(xN),max(y/sum(y)));
text(p(1),max(xN),str) % max(y/sum(y)) % max(xN) 0.15
subplot(1,2,2)
% 6&7 items
% [yN,edges] = histcounts(yr,linspace(-9,-7,31),'Normalization','probability');
% p = (edges(1:end-1)+edges(2:end))/2;

pd = fitdist(yr','Normal');
y = pdf(pd,p);
plot(p,yN,'b.',p,y/sum(y),'r--')
xlabel('Position on y axis')
ylabel('Probability')
% str = [{['Width=',num2str(pd.sigma)]},{['Height=',num2str(max(yN))]},{['FitHeight=',num2str(max(y/sum(y)))]}];
str = sprintf('Width=%.2f\nHeight=%.2f\nFitHeight=%.2f',pd.sigma,max(yN),max(y/sum(y)));
text(p(1),max(yN),str) % max(y/sum(y)) % max(yN) 0.15

%%
xrp = xr/31*pi;
yrp = yr/31*pi;
lowx = x/31*pi-pi;
upx = lowx + 2*pi;
lowy = y/31*pi-pi;
upy = lowy + 2*pi;
if xrp < lowx
    xrp = xrp + 2*pi;
elseif xrp > upx
    xrp = xrp - 2*pi;
end
if yrp < lowy
    yrp = yrp + 2*pi;
elseif yrp > upy
    yrp = yrp - 2*pi;
end
edgesx = linspace(lowx,upx,51); 
edgesy = linspace(lowy,upy,51); 
[xN,edgesx] = histcounts(xrp,edgesx,'Normalization','probability');
[yN,edgesy] = histcounts(yrp,edgesy,'Normalization','probability');
px = (edgesx(1:end-1)+edgesx(2:end))/2;
py = (edgesy(1:end-1)+edgesy(2:end))/2;
subplot(1,2,1)
pdx = fitdist(xrp','Normal');
y = pdf(pdx,px);
plot(px,xN,'b.',px,y/sum(y),'r--')
xlabel('Position on x axis')
ylabel('Probability')
str = sprintf('Width=%.2f\nHeight=%.2f\nFitHeight=%.2f',pdx.sigma,max(xN),max(y/sum(y)));
text(px(1),max(xN),str) % max(y/sum(y)) % max(xN) 0.15
subplot(1,2,2)
pdy = fitdist(yrp','Normal');
y2 = pdf(pdy,py);
plot(py,yN,'b.',py,y2/sum(y2),'r--')
xlabel('Position on y axis')
ylabel('Probability')
% str = [{['Width=',num2str(pd.sigma)]},{['Height=',num2str(max(yN))]},{['FitHeight=',num2str(max(y/sum(y)))]}];
str = sprintf('Width=%.2f\nHeight=%.2f\nFitHeight=%.2f',pdy.sigma,max(yN),max(y2/sum(y2)));
text(py(1),max(yN),str)
%%
Width = zeros(1,7);
Height = zeros(1,7);
FitHeight = zeros(1,7);
WidthWithin2s = zeros(1,7);
HeightWithin2s = zeros(1,7);
FitHeightWithin2s = zeros(1,7);
%%
item = 7;
Width(item) = (pdx.sigma + pdy.sigma)/2
Height(item) = (max(xN) + max(yN))/2
FitHeight(item) = (max(y/sum(y)) + max(y2/sum(y2)))/2
%%
item = 7;
WidthWithin2s(item) = (pdx.sigma + pdy.sigma)/2
HeightWithin2s(item) = (max(xN) + max(yN))/2
FitHeightWithin2s(item) = (max(y/sum(y)) + max(y2/sum(y2)))/2
%% Relative precision
Width = Width/Width(end);
WidthBin10ms = WidthBin10ms/WidthBin10ms(end);
WidthWithin2s = WidthWithin2s/WidthWithin2s(end);
Height = Height/Height(1);
HeightBin10ms = HeightBin10ms/HeightBin10ms(1);
HeightWithin2s = HeightWithin2s/HeightWithin2s(1);
FitHeight = FitHeight/FitHeight(1);
FitHeightBin10ms = FitHeightBin10ms/FitHeightBin10ms(1);
FitHeightWithin2s = FitHeightWithin2s/FitHeightWithin2s(1);
%%
x = 1:7;
figure
v1 = polyfit(log10(x),log10(Width),1);
y1 = 10^v1(2)*x.^v1(1);
v2 = polyfit(log10(x),log10(WidthBin10ms),1);
y2 = 10^v2(2)*x.^v2(1);
v3 = polyfit(log10(x),log10(WidthWithin2s),1);
y3 = 10^v3(2)*x.^v3(1);
plot(x,Width,'bo',x,y1,'b-',...
    x,WidthBin10ms,'g>',x,y2,'g-',...
    x,WidthWithin2s,'r*',x,y3,'r-') % 
legend('bin4ms,within7.75s','power-law fit','bin10ms,within7.75s','power-law fit','bin4ms,within2s','power-law fit')
% plot(x,Width,'o-',x,WidthBin10ms,'>--',x,WidthWithin2s,'*-.') % 
% legend('bin4ms,within7.75s','bin10ms,within7.75s','bin4ms,within2s') % 
xlabel('Loading Items')
ylabel('Width of Precision Curve')
% v1(1)
%%
x = [1,2,4,7];
figure
v1 = polyfit(log10(x),log10(FitHeightWithin2s(x)/FitHeightWithin2s(1)),1);
x1 = 1:0.1:10;
y1 = 10^v1(2)*x1.^v1(1);
plot(x,FitHeightWithin2s(x)/FitHeightWithin2s(1),'k.','MarkerSize',10)
hold on
plot(x1,y1,'b-') %  
xlabel('Loading Items')
ylabel('Precision Index')
%%
figure
v1 = polyfit(log10(x),log10(FitHeight),1);
y1 = 10^v1(2)*x.^v1(1);
v2 = polyfit(log10(x),log10(FitHeightBin10ms),1);
y2 = 10^v2(2)*x.^v2(1);
v3 = polyfit(log10(x),log10(FitHeightWithin2s),1);
y3 = 10^v3(2)*x.^v3(1);
plot(x,FitHeight,'bo',x,y1,'b-',...
    x,FitHeightBin10ms,'g>',x,y2,'g-',...
    x,FitHeightWithin2s,'r*',x,y3,'r-') % 
legend('bin4ms,within7.75s','power-law fit','bin10ms,within7.75s','power-law fit','bin4ms,within2s','power-law fit')
% plot(x,FitHeight,'o-',x,FitHeightWithin2s,'*-.')
% legend('bin4ms,within7.75s','bin4ms,within2s')
xlabel('Loading Items')
ylabel('Height of Precision Curve')
end