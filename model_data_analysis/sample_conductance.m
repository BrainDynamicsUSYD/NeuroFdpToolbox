dbstop if error
dir_strut = dir('*0_neurosamp.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
%%
V_I = -80; % mV, reversal potential
V_E = 0;
LFP_gamma = R.LFP.LFP_gamma;
ICI = [];
Amp = [];

% % set(gcf,'color','w');
% for i = 1:no
%     tempAmp = [];
%     [p1,l1] = findpeaks(LFP_gamma(i,:));
%     [p2,l2] = findpeaks(-LFP_gamma(i,:));
%     l = min(length(l1),length(l2));
%     m = 1; % index for l2
%     for j = 1:l
%         while m <= length(l2) && l1(j) >= l2(m)
%             m = m + 1;
%         end
%         if m > length(l2)
%             break
%         end
%         tempAmp = [tempAmp p1(j) + p2(m)];
%         m = m + 1;
%     end
%     tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
%     if length(tempAmp) < length(tempICI)
%         disp('Miss Match!');
%     end
%     tempAmp = tempAmp(1:length(tempICI));
%     ICI = [ICI tempICI];
%     Amp = [Amp tempAmp];
% end
figure
for i = 1% :num_files % [1 12 17 18]
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    S = load(files{i});
    disp('done.\n');
    LFPind = [1 13 14 14 11 11 11 12]; % old
    %     LFPind = [1 13 13 14 14 11 12 12]; % new
    CorrR1 = zeros(1,8);
    CorrR2 = zeros(1,8);
    for j = 1:8
        %         EPSC = (V_I - V_E)*(S.I_AMPA(j,:) + S.I_ext(j,:))./(S.V(j,:) - V_E);
        IPSC = abs((V_E - V_I)*S.I_GABA(j,:)./(S.V(j,:) - V_I));
        oneLFP = LFP.LFP{1}(LFPind(j),:);
        %         gammaLFP = LFP_gamma(LFPind(j),:);
        tempAmp = [];
        tempICI = [];
        [p1,l1] = findpeaks(oneLFP);
        [p2,l2] = findpeaks(IPSC);
        [p3,l3] = findpeaks(-IPSC);
        for t = 3:length(l1)-1
            indTop = find(l2>=l1(t) & l2<=l1(t+1));
            if ~isempty(indTop)
                indTop = indTop(1);
                indDown = find(l3<l2(indTop));
                indDown = indDown(end);
                tempAmp = [tempAmp p2(indTop)+p3(indDown)];
                tempICI = [tempICI 0.1*(l1(t+1)-l1(t))];
            end
        end
        dat = [tempAmp',tempICI'];
        n = hist3(dat,[50 50]);
        n1 = n';
        n1(size(n,1) + 1, size(n,2) + 1) = 0;
        xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
        yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
        h = pcolor(xb,yb,n1);
        set(h, 'EdgeColor', 'none');
        h.ZData = ones(size(n1)) * -max(max(n));
        colormap(hot)
        oldcmap = colormap(gray);
        colormap( flipud(oldcmap) );
        colorbar
        [r,p] = corrcoef(tempAmp',tempICI');
        if p(2) > 10^(-4)
            pvalue = ['p = ',num2str(p(2),'%.2f')];
        else
            pvalue = ['p < 10^{-4}'];
        end
        tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
        text(0.7,0.2,tp,'Units', 'Normalized','FontSize',8)
        xlabel('|IPSC(nA)|')
        ylabel('LFP IEI(ms)')
        
        %         xlim([0 100])
        %         ylim([0 6])
        
        %         dat1 = [oneLFP',EPSC'];
        %         dat2 = [oneLFP',IPSC'];
        %         n = hist3(dat1,[50 50]);
        %         n1 = n';
        %         n1(size(n,1) + 1, size(n,2) + 1) = 0;
        %         xb = linspace(min(dat1(:,1)),max(dat1(:,1)),size(n,1)+1);
        %         yb = linspace(min(dat1(:,2)),max(dat1(:,2)),size(n,1)+1);
        %         h = pcolor(xb,yb,n1);
        %         set(h, 'EdgeColor', 'none');
        %         h.ZData = ones(size(n1)) * -max(max(n));
        %         colormap(hot)
        %         oldcmap = colormap(gray);
        %         colormap( flipud(oldcmap) );
        %         colorbar
        %         [r,p] = corrcoef(oneLFP',EPSC');
        %         if p(2) > 10^(-4)
        %             pvalue = ['p = ',num2str(p(2),'%.2f')];
        %         else
        %             pvalue = ['p < 10^{-4}'];
        %         end
        %         tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
        %         text(0.7,0.2,tp,'Units', 'Normalized','FontSize',8)
        %         xlabel('LFP(a.u)')
        %         ylabel('|EPSC(nA)|')
        %         xlim([0 100])
        %         ylim([0 6])
        %         CorrR1(j) = r(2);
        %         figure
        %         n = hist3(dat2,[50 50]);
        %         n1 = n';
        %         n1(size(n,1) + 1, size(n,2) + 1) = 0;
        %         xb = linspace(min(dat2(:,1)),max(dat2(:,1)),size(n,1)+1);
        %         yb = linspace(min(dat2(:,2)),max(dat2(:,2)),size(n,1)+1);
        %         h = pcolor(xb,yb,n1);
        %         set(h, 'EdgeColor', 'none');
        %         h.ZData = ones(size(n1)) * -max(max(n));
        %         colormap(hot)
        %         oldcmap = colormap(gray);
        %         colormap( flipud(oldcmap) );
        %         colorbar
        %         [r,p] = corrcoef(oneLFP',IPSC');
        %         if p(2) > 10^(-4)
        %             pvalue = ['p = ',num2str(p(2),'%.2f')];
        %         else
        %             pvalue = ['p < 10^{-4}'];
        %         end
        %         tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
        %         text(0.7,0.2,tp,'Units', 'Normalized','FontSize',8)
        %         xlabel('LFP(a.u)')
        %         ylabel('|IPSC(nA)|')
        %         xlim([0 100])
        %         ylim([0 20])
        %         CorrR2(j) = r(2);
        %         next = input('\t Next figure?');
        %         close all
    end
    %     histogram(EPSC)
    %     [N,edges] = histcounts(EPSC);
    %     edges = (edges(1:end-1)+edges(2:end))/2;
    %     figure
    %     semilogy(edges,N)
    %     t1 = 5e4 + 1;
    %     t2 = 5e4 + 5e3;
    %     t = t1:t2;
    %     plot(t,IPSC(t),'b',t,4*EPSC(t),'r')
    %     set(gca,'XTickLabel',[0:50:500])
    %     ylabel('PSC(nA)')
    %     xlabel('Time(ms)')
    %     legend('IPSC','4*EPSC')
    %     set(gcf,'color','w');
    %     CorrR = zeros(1,size(S.I_AMPA,1));
    %     for no = 1:size(S.I_AMPA,1)
    %     g_E = -(S.I_AMPA(no,:) + S.I_ext(no,:))./(S.V(no,:) - V_E);% uS
    %     g_I = -S.I_GABA(no,:)./(S.V(no,:) - V_I);
    %     [r,~] = corrcoef(g_E,g_I);
    %     CorrR(no) = r(2);
    %     end
    
    %     g_E = -(S.I_AMPA(:) + S.I_ext(:))./(S.V(:) - V_E)*1e3;% nS
    %     g_I = -S.I_GABA(:)./(S.V(:) - V_I)*1e3;
    %     dat = [g_E,g_I];
    %     n = hist3(dat,[100 100]);
    %     n1 = n';
    %     n1(size(n,1) + 1, size(n,2) + 1) = 0;
    %     xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    %     yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
    %     h = pcolor(xb,yb,n1);
    %     set(h, 'EdgeColor', 'none');
    %     h.ZData = ones(size(n1)) * -max(max(n));
    %     colormap(hot)
    % %     xlim([0 0.05])
    % %     ylim([0 0.16])
    %     oldcmap = colormap(gray);
    %     colormap( flipud(oldcmap) );
    %     colorbar
    %     [r,p] = corrcoef(g_E,g_I);
    %     if p(2) > 10^(-4)
    %         pvalue = ['p = ',num2str(p(2),'%.2f')];
    %     else
    %         pvalue = ['p < 10^{-4}'];
    %     end
    %     tp = {['r = ',num2str(r(2),'%.2f')],pvalue};
    %     text(0.7,0.2,tp,'Units', 'Normalized'); % ,'FontSize',8)
    %     xlabel('g_E(nS)')
    %     ylabel('g_I(nS)')
end

%%
max_lag = 200; % ms 50
dt = 0.1;
% I_exc = zeros(num_files,1e5);
% I_inh = zeros(num_files,1e5);
% for i = 1:num_files
%     fprintf('Loading RYG.mat file %s...', files{i});
%     S = load(files{i});
%     I_exc(i,:) = mean(S.I_AMPA + S.I_ext);
%     I_inh(i,:) = mean(S.I_GABA);
%     disp('done.\n');
% end
% [ac1,lags1] = xcorr(mean(I_exc(3:13:end,:)),-mean(I_inh(3:13:end,:)),round(max_lag/dt),'coeff'); % 'coeff'
% [ac2,lags2] = xcorr(mean(I_exc(7:13:end,:)),-mean(I_inh(7:13:end,:)),round(max_lag/dt),'coeff');
% [ac3,lags3] = xcorr(mean(I_exc(13:13:end,:)),-mean(I_inh(13:13:end,:)),round(max_lag/dt),'coeff');
% plot(lags1*dt,ac1,'g',lags2*dt,ac2,'r',lags3*dt,ac3,'b')
acCell = cell(1,num_files);
lagsCell = cell(1,num_files);
for id_out = 1:num_files
    fprintf('Loading RYG.mat file %s...\n', files{id_out});
    S = load(files{id_out});
    I_exc = S.I_AMPA + S.I_ext;
    I_inh = S.I_GABA;
    AC = [];
    LAGS = [];
    for i = 1:8
        [ac,lags] = xcorr(I_exc(i,:),-I_inh(i,:),round(max_lag/dt),'coeff');
        AC = [AC;ac];
        LAGS = [LAGS;lags];
    end
    acCell{id_out} = mean(AC);
    lagsCell{id_out} = mean(LAGS);
end
%%
% plot(lagsCell{3}*dt,acCell{3},'g',lagsCell{7}*dt,acCell{7},'r',lagsCell{13}*dt,acCell{13},'b')
lagsMat = cell2mat(lagsCell');
acMat = cell2mat(acCell');
% plot(mean(lagsMat(3:13:end,:))*dt,mean(acMat(3:13:end,:)),'g',mean(lagsMat(7:13:end,:))*dt,mean(acMat(7:13:end,:)),'r',mean(lagsMat(13:13:end,:))*dt,mean(acMat(13:13:end,:)),'b')
x = mean(lagsMat(3:13:end,:))*dt;
y = acMat(3:13:end,:);
shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
hold on;
x = mean(lagsMat(7:13:end,:))*dt;
y = acMat(7:13:end,:);
shadedErrorBar(x,mean(y,1),std(y),'lineprops','r');
hold on;
x = mean(lagsMat(13:13:end,:))*dt;
y = acMat(13:13:end,:);
shadedErrorBar(x,mean(y,1),std(y),'lineprops','b');
legend('IE ratio = 0.8','IE ratio = 1.0','IE ratio = 1.3')
xlabel('Lags(ms)')
ylabel('E-I current xcorr')
%%
subplot(3,2,1)
t = 5e3+(1:1e4);
plot(t*dt,I_exc1(1,t),t*dt,I_inh1(1,t))
title('IE ratio = 0.8')
xlabel('Time(ms)')
ylabel('Current(nA)')
subplot(3,2,3)
plot(t*dt,I_exc2(1,t),t*dt,I_inh2(1,t))
title('IE ratio = 1.0')
xlabel('Time(ms)')
ylabel('Current(nA)')
% legend('Excitatory','Inhibitory')
subplot(3,2,5)
plot(t*dt,I_exc3(1,t),t*dt,I_inh3(1,t))
title('IE ratio = 1.3')
xlabel('Time(ms)')
ylabel('Current(nA)')
%%
% figure(1)
% figure_width = 8.4; %cm
% figure_hight = 15; %cm
% figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', ...
%     'color','w', 'position', [0, 0, figure_width, figure_hight], ...
%     'PaperSize', [figure_width, figure_hight]); % this is the trick!
% figure('NumberTitle','off','name', 'figure_size_control', 'units', 'centimeters', ...
%     'color','w'); % this is the trick!

[N,edges] = histcounts(I_exc1(5,:),100);
Y = (edges(1:end-1)+edges(2:end))/2;
pd = fitdist(I_exc1(5,:)','Stable')
y = pdf(pd,Y);
plot(Y,N/sum(N),'LineWidth',2)
hold on
plot(Y,y/sum(y),'m-.','LineWidth',2);
xlabel('I_{exc}(nA)') %, 'fontsize', 30 )
ylabel('Probability') % , 'fontsize', 30 )
% title('IE ratio = 0.8')
axes('Position',[.5 .5 .4 .4])
box on
loglog(Y,N/sum(N),'LineWidth',2)
hold on
loglog(Y,y/sum(y),'m-.','LineWidth',2);
% xlim(minmax(I_exc1(:)'))
xlim([0.5 max(I_exc2(5,:)')])
set(gca,'TickLength',[0.05, 0.01])
set(gca,'XTickLabel',[],'YTickLabel',[])

% text(2.5,2,{'The figure size should',' be 8.4 by 15 cm'})
% set(gca,'FontSize', 20);  % 6 points for x-axis tickmark labels
% xlabel('30 point label', 'fontsize', 30 ); % this must be after the above line!

% set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
% print -depsc figure_size_control % this is the trick!!

% subplot(3,2,2);
%%
% h2 = subplot(3,2,4);
figure(2)
[N,edges] = histcounts(I_exc2(4,:),100);
Y = (edges(1:end-1)+edges(2:end))/2;
pd = fitdist(I_exc2(4,:)','Stable')
y = pdf(pd,Y);
plot(Y,N/sum(N),'LineWidth',2)
hold on
plot(Y,y/sum(y),'m-.','LineWidth',2);
xlabel('I_{exc}(nA)')
ylabel('Probability')
% title('IE ratio = 1.0')
axes('Position',[.5 .5 .4 .4])
box on
loglog(Y,N/sum(N),'LineWidth',2)
hold on
loglog(Y,y/sum(y),'m-.','LineWidth',2);
xlim([0.5 max(I_exc2(4,:))])
set(gca,'TickLength',[0.05, 0.01])
set(gca,'XTickLabel',[],'YTickLabel',[])
%%
% h3 = subplot(3,2,6);
figure(3)
[N,edges] = histcounts(I_exc3(1,:),100);
Y = (edges(1:end-1)+edges(2:end))/2;
pd = fitdist(I_exc3(1,:)','Stable')
y = pdf(pd,Y);
plot(Y,N/sum(N),'LineWidth',2)
hold on
plot(Y,y/sum(y),'m-.','LineWidth',2);
xlabel('I_{exc}(nA)')
ylabel('Probability')
% title('IE ratio = 1.3')
axes('Position',[.5 .5 .4 .4])
box on
loglog(Y,N/sum(N),'LineWidth',2)
hold on
loglog(Y,y/sum(y),'m-.','LineWidth',2);
xlim([0.5 max(I_exc3(1,:))])
set(gca,'TickLength',[0.05, 0.01])
set(gca,'XTickLabel',[],'YTickLabel',[])