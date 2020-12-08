% function DistributionDisplacement
% displacement distribution
% adapt from Xian's MSDnew.m function
% close all
dir_strut = dir('3DBurst3000*minTime30SR1000.mat'); % 3DBurstLFP0*minTime30SR1000P95.mat
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
hw = 31;
Dis = [];
deltaTime = [2 4 8 16 32]; % 30 ;
figure
% subplot(2,2,[3 4])
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
%             posCenter = posCenter/hw*pi;
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            Displace = [Displace;DisplaceTemp] ;
        end
    end
    %     rmsd = sqrt(mean(Displace));
%     Displace = sqrt(Displace);
%     rmsd = deltaT^eta;
    Displace = Displace/rmsd;
    [n,x] = histcounts(Displace,50,'normalization','probability') ; % 80 50
%     plot(x(2:end),n/sum(n),'.','MarkerSize',6)
%     loglog((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    semilogy((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    hold on;
    Dis = [Dis; Displace];
%     if deltaT == 8 % 16
%         %         pd = fitdist(Displace,'normal'); % normal
%         %         y = pdf(pd,x) ;
%         x1 = x;
%     end
%     if deltaT == 8 % 8
% %                 pd = fitdist(Displace,'Stable')
% %                 yy = pdf(pd,x);
%         x2 = x;
%     end
    %     if deltaT == 32
    %         loglog(x1,y/sum(y),'r--','LineWidth',2); % 0.1*
    %         hold on ;
    %         loglog(x2,yy/sum(yy),'m-.','LineWidth',2);
    %     end
    
    %     legend(['t = ', num2str(deltaT),' ms'],'Gaussian')
    
    %     pd = fitdist(Displace,'HalfNormal')
    %     y = pdf(pd,x) ;
    %     hold on ; loglog(x,10^3.5*y,'r--');
    %         legend('t = 30 ms','Gaussian')
    %     xlabel('r(t)')
    %     ylabel('P(r(t))')
    %     next = input('\t Next figure?');
    %     delete(gcf);
end
% legend({'\tau = 2ms','\tau = 4 ms','\tau = 8 ms','\tau= 16 ms','\tau= 32 ms','Gaussian','\alpha stable'},'fontsize',8)
xlabel('Displacement d(\tau)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d(\tau))','fontsize',8)
ax=axes('Position',[0.65 0.65 0.2 0.2],'Unit','normalize',...
    'parent',1);
box on
%%
hw = 31;
Dis = [];
eta = 0.6; % 0.9 % 0.67
deltaTime = [2 4 8 16 32]; % 30 ;
figure
% subplot(2,2,[3 4])
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
%             posCenter = posCenter/hw*pi;
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            Displace = [Displace;DisplaceTemp] ;
        end
    end
    %     rmsd = sqrt(mean(Displace));
%     Displace = sqrt(Displace);
    rmsd = deltaT^eta;
    Displace = Displace/rmsd;
    [n,x] = histcounts(Displace,50,'normalization','probability') ; % 80 50
%     plot(x(2:end),n/sum(n),'.','MarkerSize',6)
%     loglog((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    semilogy((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    hold on;
    Dis = [Dis; Displace];
%     if deltaT == 8 % 16
%         %         pd = fitdist(Displace,'normal'); % normal
%         %         y = pdf(pd,x) ;
%         x1 = x;
%     end
%     if deltaT == 8 % 8
% %                 pd = fitdist(Displace,'Stable')
% %                 yy = pdf(pd,x);
%         x2 = x;
%     end
    %     if deltaT == 32
    %         loglog(x1,y/sum(y),'r--','LineWidth',2); % 0.1*
    %         hold on ;
    %         loglog(x2,yy/sum(yy),'m-.','LineWidth',2);
    %     end
    
    %     legend(['t = ', num2str(deltaT),' ms'],'Gaussian')
    
    %     pd = fitdist(Displace,'HalfNormal')
    %     y = pdf(pd,x) ;
    %     hold on ; loglog(x,10^3.5*y,'r--');
    %         legend('t = 30 ms','Gaussian')
    %     xlabel('r(t)')
    %     ylabel('P(r(t))')
    %     next = input('\t Next figure?');
    %     delete(gcf);
end
% x1 = [linspace(-60,-50,100) x1 linspace(50,60,100)];
% x2 = [linspace(-200,-80,100) x2 linspace(80,200,100)];
Dis = [-Dis;Dis];

hold on;
pd = fitdist(Dis,'normal') % normal
% x1 = linspace(0,0.35*max(Dis),130);
x1 = linspace(-0.35*max(Dis),0.35*max(Dis),130);
y = pdf(pd,x1); % x1
pd = fitdist(Dis,'Stable')
x2 = linspace(-max(Dis),max(Dis),130);
yy = pdf(pd,x2);
% loglog(x1,y/sum(y),'k--','LineWidth',2); % 0.1*
semilogy(x1,y/sum(y),'k--','LineWidth',2); %*(x1(2)-x1(1))
hold on ;
% loglog(x2,yy/sum(yy),'r-','LineWidth',2);
semilogy(x2,yy/sum(yy),'r-','LineWidth',2); % *(x2(2)-x2(1))
% xlim([0.9 1e3])

legend({'\tau = 2ms','\tau = 4 ms','\tau = 8 ms','\tau= 16 ms','\tau= 32 ms','Gaussian','\alpha stable'},'fontsize',8,'edgecolor','none','color','none')
xlabel('Displacement d_s(\tau)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d_s(\tau))','fontsize',8)
% xlabel('Displacement d(\tau)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
% ylabel('{\itP} (d(\tau))','fontsize',8)

% xlabel(['displacement/\eta^{', num2str(eta),'}']) % ('displacement/\eta^{0.5}')
% ylabel('Probability')
% title('Displacement Distribution')
% text(-0.08,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% rmsd = sqrt(mean(Displace));
% Displace = sqrt(Displace);
% % Displace = Displace(Displace > 80);
% Displace = Displace/rmsd;
% [n,x] = hist(Displace,240) ;
% loglog(x,n/sum(n),'o')
% hold on
% legendInfo{delta} = ['t = ', num2str(deltaT),' ms'] ;
% pd = fitdist(Displace,'HalfNormal') % normal
% y = pdf(pd,x) ;
% hold on ; loglog(x,y,'r--'); % 1.2e3*
% legend('t = 30 ms','Gaussian')
% xlabel('r(t)')
% ylabel('P(r(t))')
% % end
%
ax=axes('Position',[0.65 0.65 0.2 0.2],'Unit','normalize',...
    'parent',1);
box on
Dis = [];
deltaTime = [2 4 8 16 32]; % 30 ;
for delta = 1:length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            Displace = [Displace;DisplaceTemp] ;
        end
    end
    %     rmsd = sqrt(mean(Displace));
%     Displace = sqrt(Displace);
%     rmsd = deltaT^eta;
    Displace = Displace/rmsd;
    [n,x] = histcounts(Displace,50,'normalization','probability') ; % 80 50
%     plot(x(2:end),n/sum(n),'.','MarkerSize',6)
%     loglog((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    semilogy((x(1:end-1)+x(2:end))/2,n/sum(n),'.','MarkerSize',6)
    hold on;
    Dis = [Dis; Displace];
end
xlabel('Displacement d(\tau)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d(\tau))','fontsize',8)
%% adjust spikes pattern center data type
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
minBurstTime = 0;
fs = 1e3;
for i = 1:num_files % (num_files - 12):num_files
    %     loop_num = loop_num + 1;
    %     % For PBS array job
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    raw = ~isnan(R.grid.bayes.radius); % .grid
    ind = find(R.grid.bayes.bayes_factor_ln <= log(100)); % .grid
    raw(ind) = 0;
    t_mid = R.grid.t_mid;
    centre = R.grid.bayes.centre;
    [~,high_du1,~,start,~] = seq_postprocess(raw,1); % ms
    WCentroids = cell(1,length(high_du1));
    for j = 1:length(high_du1)
        WCentroids{j} = [t_mid(start(j):(start(j)+high_du1(j)-1))' centre(:,start(j):(start(j)+high_du1(j)-1))'];
    end
    % 2DBurst4: spikes pattern based spatial pattern
    saveFileName = ['3DBurst4',sprintf('%04g',i),...
        'minTime',num2str(minBurstTime),'SR',num2str(fs),'.mat'] ;
    save(saveFileName,'WCentroids') ;
end
%% generate random walk
% one event [t,x,y]
% R.WCentroids = cell(1);
tot = 1e4;
R.WCentroids{1} = zeros(tot,3);
R.WCentroids{1}(:,1) = 1:tot;
x = 0;
y = 0;
for i = 1:tot
    r = rand;
    alpha = rand*2*pi;
    x = x + rand*cos(alpha);
    y = y + rand*sin(alpha);
    R.WCentroids{1}(i,2:3) = [x,y];
end
%% multiple events
t0 = round(100*rand);
for j = 1:200
    tot = 30+round(30*rand);
    R.WCentroids{j} = zeros(tot,3);
    t0 = t0 + round(120*rand);
    R.WCentroids{j}(:,1) = t0+(1:tot);
    x = 0;
    y = 0;
    for i = 1:tot
        r = rand;
        alpha = rand*2*pi;
        x = x + rand*cos(alpha);
        y = y + rand*sin(alpha);
        R.WCentroids{j}(i,2:3) = [x,y];
    end
end
%%
hw = 31;
deltaTime = [2 4 8 16 32]; % 30 ;
Color = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880];
subplot(4,2,[7,8])
for delta = 1:3 % length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            DisplaceTemp = sqrt(sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2));
%             DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            Displace = [Displace;DisplaceTemp] ;
        end
    end
%     nEdge = logspace(log10(min(Displace(Displace>0))),log10(max(Displace)),20) ;
%     [n,nEdge] = histcounts(Displace,nEdge,'normalization','pdf') ;
%     x = 10.^((log10(nEdge(1:end-1))+log10(nEdge(2:end)))/2);
%     [n,nEdge] = histcounts(log10(Displace),40,'normalization','pdf');
%     x = 10.^((nEdge(1:end-1)+nEdge(2:end))/2);
    [n,nEdge] = histcounts(Displace,80,'normalization','probability');
    x = (nEdge(1:end-1)+nEdge(2:end))/2;
    h(delta) = loglog(x,n/sum(n),'.','MarkerSize',6,'color',Color(delta,:));
    hold on;
    body = Displace(Displace<200);
    pd = fitdist([-body;body],'Stable')
    xbody = x(x<200);
    y = pdf(pd,xbody);
    loglog(xbody,y/sum(y),'-','LineWidth',1,'color',Color(delta,:));
    hold on ;
end
xlabel('Displacement d(\tau) (\mum)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d(\tau))','fontsize',8)
legend([h(1) h(2) h(3)],{'\tau = 2 ms','\tau = 4 ms','\tau = 8 ms',},'fontsize',8,'edgecolor','none','color','none')
text(-0.2,1.02,'E','Units', 'Normalized','FontSize',12)
ax=axes('Position',[0.85 0.35 0.1 0.1],'Unit','normalize',...
    'parent',1);
box on

Dis = [];
eta = 0.6; % 0.9 % 0.67
for delta = 1:3 %length(deltaTime)
    deltaT = deltaTime(delta) ;
    Displace = [] ;
    for i = 1:num_files % 7:13:num_files % 1:num_files
        fprintf('Loading 3DBurst.mat file %s...\n', files{i});
        R = load(files{i});
        center = R.WCentroids ;
        for iBurst = 1:size(center,2)
            posCenter = center{iBurst}(:,2:3);
            posCenter = posCenter/40*600; %% manually modify here according to electrdoes %%
            posCenter = posCenter/600*2*pi;
            posCenter = exp(1i*posCenter);
%             DisplaceTemp = angle(posCenter(deltaT+1:end,1)./posCenter(1:end-deltaT,1))/(2*pi)*600;
            DisplaceTemp = sqrt(sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2));
%             DisplaceTemp = sum((angle(posCenter(deltaT+1:end,:)./posCenter(1:end-deltaT,:))/(2*pi)*600).^2,2);
            Displace = [Displace;DisplaceTemp] ;
        end
    end
    rmsd = deltaT^eta;
    Displace = Displace/rmsd;
    [n,nEdge] = histcounts(Displace,80,'normalization','probability');
    x = (nEdge(1:end-1)+nEdge(2:end))/2;
    if delta == 1
        x1 = x;
    end
%     if delta == 2
%         Dis2 = Displace;
%     end
    loglog(x,n/sum(n),'.','MarkerSize',6,'color',Color(delta,:))        
    hold on ;
    Dis = [Dis; Displace];
end
body = Dis;
pd = fitdist([-body;body],'Stable')
xbody = [x,max(x)+1:max(x1)];
y = pdf(pd,xbody);
h(1) = loglog(xbody,y/sum(y),'r-','LineWidth',1);
pd = fitdist([-body;body],'normal')
xbody = x(x<100);
y = pdf(pd,xbody);
h(2) = loglog(xbody,y/sum(y),'k--','LineWidth',1);
xlabel('Displacement d(\tau) (\mum)','fontsize',8) % ('displacement/\eta^{0.5}') ['displacement/\eta^{', num2str(eta),'}']
ylabel('{\itP} (d(\tau))','fontsize',8)
legend([h(1) h(2)],{'\alpha stable','Gaussian'},'fontsize',8,'edgecolor','none','color','none')