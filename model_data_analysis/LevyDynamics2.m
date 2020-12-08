
% plot the fitted levy distribution of jump
% cd /scratch/RDS-FSC-cortical-RW/longtrialforlevy
clear
close all
d = dir('*RYG.mat');
%g_balance = 0.7:0.01:1.3;
repeats = 10;
for a=50 %0:10:200
    window_T = a; % unit : 0.1 ms
    jump_win= a;
    for ii = 13 % :length(d)/repeats
        clear distance_spon delta_x_spon
        kk=1;
        for jj = ii:length(d)/repeats:length(d) % (ii-1)*repeats+1:ii*repeats
            R=load(d(jj).name); %,'grid','neuron_stats');
            %             R = get_grid_firing_centreYL(R);
            [distance_spon{kk},delta_x_spon{kk},factor{kk}] = FitLevyDynamics_fun2(R,'bayes',window_T,jump_win); % bayes
            %         [distance_spon{kk},delta_x_spon{kk}] = FitLevyDynamics_fun(R,'real-time_I');
            kk=kk+1;
        end
        close all
        % edge_custom =linspace(0,max(cellfun(@max,distance_spon)),21);
        dis = cell2mat(distance_spon');
        A = cell2mat(delta_x_spon');
%                 factor = cell2mat(factor);
%                 factor(factor <= log(100)) = 0;
%                 factor(factor > 0) = 1;
%                 ff = factor(1:end-jump_win/10).*factor((jump_win/10+1):end);
%                 A = A(ff>0);
%                 dis = dis(ff>0);
        
        %         A = A(~isnan(A));
        %         dis = dis(~isnan(dis));
        f = fitdist([A(:);-A(:)],'stable')
        %         f2 = fitdist([dis(:)],'stable');
        f2 = fitdist([dis(:);-dis(:)],'stable')
        hw=31;%pi
        for  bin_num =30%11:5:101
            edge_custom =linspace(0,sqrt(2)*hw,bin_num);
            clear N
            for ll =1:repeats
                [N(:,ll),edge]=histcounts(distance_spon{ll},51,'normalization','pdf');
            end
            figure('visible','on')
            subplot(2,2,1)
            shadedErrorBar(edge(2:end),mean(N,2),std(N,[],2))
            title('Distance')
            xlabel('Distance (grid points)')
            ylabel('Probability')
            set(gca,'yscale','log','xscale','log')
            set(gca,'linewidth',1.5,'fontsize',15)
            
            subplot(2,2,2)
            [N,edge]=histcounts(dis(:),51,'normalization','pdf');
            loglog(edge(2:end),N,'.')
            title('Distance')
            xlabel('Distance (grid points)')
            ylabel('Probability')
            set(gca,'linewidth',1.5,'fontsize',15)
            % histogram(distance,'normalization','probability');
            
            
            subplot(2,2,3)
            cla
            x = linspace(-hw,hw,1000);
            plot(x,pdf(f,x),'r-','linewidth',2)
            hold on
            
            % plot(edge(1:end-1),N)
            h=histogram([A(:);-A(:)],'normalization','pdf');% #bins:60
            xlabel('X change (grid points)')
            ylabel('Probability')
            %         legend({'fitted','data'})
            set(gca,'linewidth',1.5,'fontsize',15)
            
            subplot(2,2,4)
%                         cla
%                         x = linspace(-hw,hw,1000);
%                         plot(x,pdf(f2,x),'r-','linewidth',2)
%                         hold on
%                         h=histogram([dis(:);-dis(:)],'normalization','pdf');% #bins:60
%                         xlabel('Distance (grid points)')
%                         ylabel('Probability')
%                         set(gca,'linewidth',1.5,'fontsize',15)
            
            logbin = logspace(-2,0.5,bin_num);%dlogbin = diff(logbin);
            h=histogram(dis(:),logbin,'normalization','pdf');% #bins:60
            set(gca,'linewidth',1.5,'fontsize',15)
            %saveas(gcf,sprintf('g_%.2f_win_%d.jpg',g_balance(ii),a))
%             saveas(gcf,sprintf('pop_levy_dist_tau_%.1f_win_%d.fig',ii,a/10))
%             save(sprintf('pop_levy_dist_tau_%.1f_win_%d.mat',ii,a/10),'A','dis','f','-v7.3')
        end
    end
end

% for a = 10:10:200
% window_T = a; % unit : 0.1 ms
% jump_win= a;
% % fit in stable distribution
% [distance_spon,delta_x_spon,distance_evok,delta_x_evok] = deal(cell(1,10));%length(d)));
% jj = 1;
% for ii=1:10%length(d)
%     R=load(d(ii).name);
%     [distance_spon{jj},delta_x_spon{jj}] = FitLevyDynamics_fun(R,'pop',window_T,jump_win);
% %     [pd_evok{ii},distance_evok{ii},delta_x_evok{ii},x_variable_evok{ii},y_variable_evok{ii}] = FitLevyDynamics_fun(R,'populationvector',floor(4e4/window_T):floor(1e5/window_T-20),0,window_T/10);
%     jj=jj+1;
% end
% save(sprintf('LevyDynamics_%d.mat',window_T),'-v7.3')
% % load('LevyDynamics_20.mat')
% % plotLevyDynamics(pd_spon,distance_spon,delta_x_spon,x_variable_spon,y_variable_spon,'spon',bins)
% % plotLevyDynamics(pd_evok,distance_evok,delta_x_evok,x_variable_evok,y_variable_evok,'evok',bins)
%
%
% close all
% % edge_custom =linspace(0,max(cellfun(@max,distance_spon)),21);
% dis = cell2mat(distance_spon');
% A = cell2mat(delta_x_spon');
% f = fitdist(A(:),'stable');
% hw=31;%pi
% for  bin_num =51%11:5:101
%     edge_custom =linspace(0,sqrt(2)*hw,bin_num);
%     clear N
%     for ii =1:10
%         [N(:,ii),edge]=histcounts(distance_spon{ii},51,'normalization','pdf');
%     end
%     figure('visible','off')
%     subplot(2,2,1)
%     shadedErrorBar(edge(2:end),mean(N,2),std(N,[],2))
%     title('Distance')
%     xlabel('Distance (grid points)')
%     ylabel('Probability')
%     set(gca,'yscale','log','xscale','log')
%
%     subplot(2,2,2)
%     [N,edge]=histcounts(dis(:),51,'normalization','pdf');
%     loglog(edge(2:end),N,'.')
%     title('Distance')
%     xlabel('Distance (grid points)')
%     ylabel('Probability')
%
%     % histogram(distance,'normalization','probability');
%
%
%     subplot(2,2,3)
%     x = linspace(-hw,hw,bin_num);
%     plot(x,pdf(f,x),'r-','linewidth',2)
%     hold on
%
%     % plot(edge(1:end-1),N)
%     edge_custom =linspace(-hw,hw,bin_num);
%     h=histogram(A,51,'normalization','pdf');% #bins:60
%     xlabel('X change (grid points)')
%     ylabel('Probability')
%     legend({'fitted','data'})
%     title(sprintf('X change,bin:%d',bin_num))
%
%     subplot(2,2,4)
%     logbin = logspace(-2,0.5,bin_num);%dlogbin = diff(logbin);
%     h=histogram(dis(:),logbin,'normalization','pdf');% #bins:60
%     saveas(gcf,sprintf('bin_%d.jpg',a))
% end
% end
%
% function plotLevyDynamics(pd,distance,delta_x,x_variable,y_variable,state,bins)
% g_balance = (0.7:0.05:1.3);
% % plot
% figure_width = 16;
% figure_hight = 10.88;
% total_row = 2;
% total_column = 3;
% SV = 0.1176;
% SH = 0.08;
% figure('NumberTitle','off','name', 'Reproduction', 'units', 'centimeters', ...
%     'color','w', 'position', [0, 0, figure_width, figure_hight], ...
%     'PaperSize', [figure_width, figure_hight]); % this is the trick!
% offset =1;
% % 0.85
% subaxis(total_row,total_column,1,1,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% plot(x_variable{(4-1)*10+offset},y_variable{(4-1)*10+offset},'r-','LineWidth',1.5)
% hold on
% if nargin > 6
%     h=histogram(delta_x{(4-1)*10+offset},bins,'normalization','probability');
% else
%     h=histogram(delta_x{(4-1)*10+offset},'normalization','probability');
% end
% % axis([-30 30 -inf inf])
% xlabel('\Delta x (grid point)','FontSize',10,'FontWeight','bold');
% y=ylabel('p(\Delta x)','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% title('g/g_c=0.85','FontSize',10,'FontWeight','bold');
% text(-0.25,1.15,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
%
% % 1
% subaxis(total_row,total_column,2,1,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% plot(x_variable{(7-1)*10+offset},y_variable{(7-1)*10+offset},'r-','LineWidth',1.5)
% hold on
% h=histogram(delta_x{(7-1)*10+offset},bins,'normalization','probability');
% % axis([-30 30 -inf inf])
% xlabel('\Delta x (grid point)','FontSize',10,'FontWeight','bold');
% y=ylabel('p(\Delta x)','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% title('g/g_c=1','FontSize',10,'FontWeight','bold');
% % text(-0.25,1.15,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
%
% % 1.15
% subaxis(total_row,total_column,3,1,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% plot(x_variable{(10-1)*10+offset},y_variable{(10-1)*10+offset},'r-','LineWidth',1.5)
% hold on
% h=histogram(delta_x{(10-1)*10+offset},bins,'normalization','probability');
% % axis([-30 30 -inf inf])
% xlabel('\Delta x (grid point)','FontSize',10,'FontWeight','bold');
% y=ylabel('p(\Delta x)','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% title('g/g_c=1.15','FontSize',10,'FontWeight','bold');
% % text(-0.25,1.15,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
%
% % distance
% subaxis(total_row,total_column,1,2,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% [N,edge]=histcounts(distance{(7-1)*10+offset},'normalization','probability');
% loglog(edge(1:end-1),N,'k-','LineWidth',1.5)
% hold on
% xrange = logspace(0.01,4,100); %xvalues for added power law slope
% alpha = -2.44; %mean field exponent...1.5
% g= 10^-0.1*(xrange.*10^0) .^ (alpha); % y values
% plot(xrange,g,'--r','linewidth',1.5);
% % axis([0.5 100 10^-5 1])
% xlabel('|\Delta r| (grid point)','FontSize',10,'FontWeight','bold');
% y=ylabel('p(|\Delta r|)','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% % title('g/g_c=0.85','FontSize',10,'FontWeight','bold');
% text(-0.25,1.15,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
%
%
% % get parameters and their confidence interval
% for ii = 1:length(pd)/10
%     %     CI = paramci(pd{ii});
%     %     alpha_CI(:,ii) =CI(:,1);
%     %     beta_CI(:,ii) = CI(:,2);
%     %     gama_CI(:,ii) = CI(:,3);
%     %     delta_CI(:,ii) = CI(:,4);
%     kk=1;
%     for jj = ((ii-1)*10+1):((ii)*10)
%         alpha(kk,ii) = pd{jj}.alpha;
%         beta(kk,ii) = pd{jj}.beta;
%         gama(kk,ii) = pd{jj}.gam;
%         delta(kk,ii) = pd{jj}.delta;
%         kk=kk+1;
%     end
%     % alpha_CI = abs(flipud(alpha_CI)-alpha);
%     % beta_CI = abs(flipud(beta_CI)-beta);
%     % gama_CI = abs(flipud(gama_CI)-gama);
%     % delta_CI = abs(flipud(delta_CI)-delta);
% end
% % alpha
% subaxis(total_row,total_column,2,2,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% shadedErrorBar(g_balance,alpha,{@mean,@std},'-k',1);
% % axis([-10 30 -10 30])
% xlabel('g/g_c','FontSize',10,'FontWeight','bold');
% y=ylabel('\alpha','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% % title('g/g_c=0.85','FontSize',10,'FontWeight','bold');
% text(-0.25,1.15,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
%
% % beta
% subaxis(total_row,total_column,3,2,'SpacingHoriz',SH,...
%     'SpacingVert',SV,'MR',0.5*SH,'ML',SH,'MT',0.5*SV,'MB',SV);
% shadedErrorBar(g_balance,beta,{@mean,@std},'-k',1);
% % axis([-10 30 -10 30])
% xlabel('g/g_c','FontSize',10,'FontWeight','bold');
% y=ylabel('\beta','FontSize',10,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% ytickangle(45)
% % title('g/g_c=0.85','FontSize',10,'FontWeight','bold');
% % text(-0.25,1.15,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',14,'FontWeight','bold')
% set(gcf, 'PaperPositionMode', 'auto');
% cd /scratch/RDS-FSC-cortical-RW/outputimage
% print(['LevyDyanmics',state,'.eps'],'-depsc','-r300')
% print(['LevyDyanmics',state,'.tif'],'-dtiff','-r300')
% % figure
% % subplot(1,2,1)
% % shadedErrorBar(g_balance,gama,{@mean,@std},'-k',1);
% % xlabel('g/g_c','FontSize',10,'FontWeight','bold');
% % ylabel('\gama','FontSize',10,'FontWeight','bold');
% % set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% %
% % subplot(1,2,2)
% % shadedErrorBar(g_balance,delta,{@mean,@std},'-k',1);
% % xlabel('g/g_c','FontSize',10,'FontWeight','bold');
% % ylabel('\delta','FontSize',10,'FontWeight','bold');
% % set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% % set(gcf, 'PaperPositionMode', 'auto');
% % print('LevyDyanmics.eps','-depsc','-r300')
% % print('LevyDyanmics.tif','-dtiff','-r300')
%
% % distance distribution
% figure
% h=histogram(distance{(7-1)*10+1},'normalization','probability');
% % axis([-10 30 -10 30])
% xlabel('|\Delta r| (grid point)','FontSize',10,'FontWeight','bold');
% ylabel('p(|\Delta r|)','FontSize',10,'FontWeight','bold');
% title(state,'FontSize',10,'FontWeight','bold');
% set(gca,'FontSize',10,'Fontname', 'Arial','FontWeight','bold','LineWidth',1.5)
% set(gcf, 'PaperPositionMode', 'auto');
% cd /scratch/RDS-FSC-cortical-RW/outputimage
% print(['Distance',state,'.eps'],'-depsc','-r300')
% print(['Distance',state,'.tif'],'-dtiff','-r300')
% end












