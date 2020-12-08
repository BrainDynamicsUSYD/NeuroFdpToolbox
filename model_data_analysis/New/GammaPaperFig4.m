% Gamma Pattern Paper Figure 4
figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'GammaPaperFig4V3', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
FR = zeros(1,130); % num_files);
for i = 1:130 % num_files
    fprintf('Processing output file No.%d out of %d...\n', i, num_files);
    fprintf('\t File name: %s\n', files{i});
    R = load(files{i});
    FR(i) = sum(R.num_spikes{1})/3969/10; % 20s or 10s
end
FR = vec2mat(FR,13);
fr = mean(FR);
frSTD = std(FR);
Fr = 3.375*(0.7:0.05:1.3);
v1 = polyfit(log10(Fr(1:7)),log10(fr(1:7)),1);
x1 = Fr(1:9);
y1 = 10^v1(2)*x1.^v1(1);
v2 = polyfit(log10(Fr(7:end)),log10(fr(7:end)),1);
x2 = Fr(5:end);
y2 = 10^v2(2)*x2.^v2(1);

subplot(2,2,1)
errorbar(Fr(1:6),fr(1:6),frSTD(1:6),'rs','MarkerSize',8,'CapSize',8,'LineWidth',1.5)
hold on
errorbar(Fr(7),fr(7),frSTD(7),'k>','MarkerSize',8,'CapSize',8,'LineWidth',1.5)
hold on
errorbar(Fr(end-5:end),fr(end-5:end),frSTD(end-5:end),'bo','MarkerSize',8,'CapSize',8,'LineWidth',1.5)
set(gca,'YScale','log')
hold on
semilogy(x1,y1,'r','LineWidth',2)
hold on
semilogy(x2,y2,'b','LineWidth',2)
hold on
x = [0.7 0.95 1.05 1.3 1.3 1.05 0.95 0.7]*3.375;
y = [ones(1,4) 100*ones(1,4)];
v= [x' y'];
f = [1 2 7 8;2 3 6 7;3:6];
col = [0.3;1;0.7];
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.3)
xlim([0.7 1.3]*3.375)
str = {'I','II','III'};
text([0.82 0.99 1.16]*3.375,[50 27 5],str,'FontSize',12)
text(-0.28,1,'A','Units', 'Normalized','FontSize',12)
xlabel('I-E ratio \xi','FontSize',10)
ylabel('Firing Rate (Hz)','FontSize',10)
% axes('Position',[.5 .5 .35 .35])
% box on
% errorbar(Fr,fr,frSTD,'o') % ,'MarkerSize',8,'CapSize',8,'LineWidth',1.5)
% hold on
% plot(x1,y1) %,'LineWidth',2)
% hold on
% plot(x2,y2) % ,'LineWidth',2)
% hold on
% patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.3)
% xlim([0.7 1.3])

hw = 31;
fw = 2*hw+1;
mode = 'quick';
[Lattice, ~] = lattice_nD(2, hw);
ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];
id = [7 3 13];
str2 = {'B','C','D'};
str3 = {'I-E ratio \xi=3.4','I-E ratio \xi=2.7','I-E ratio \xi=4.4'};
tp = cell(1,3);

tp{1} = 3750:3902;
tp{2} = 200:352;
tp{3} = 650:802;
for id_out = 2:4
    subplot(2,2,id_out)
    R = load(files{id(id_out-1)});
    t_mid = R.grid.t_mid;
    ind_a_vec = R.grid.ind_ab(1,:);
    ind_b_vec = R.grid.ind_ab(2,:);
    switch mode
        case 'bayesian'
            x_centre = R.grid.bayes.centre(1,:);
            y_centre = R.grid.bayes.centre(2,:);
            width = R.grid.bayes.radius;
        case 'quick'
            x_centre = R.grid.quick.centre(1,:);
            y_centre = R.grid.quick.centre(2,:);
            width = R.grid.quick.radius;
    end
    x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
    y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
    axis equal;
    box on;
%     set(gca,'xtick',[],'ytick',[]);
    xlabel('X (\mum)')
    ylabel('Y (\mum)')
    title(str3{id_out-1})
    xlim([-hw hw]);
    ylim([-hw hw]);
    set(gca,'xtick',[-hw 0 hw],'ytick',[-hw 0 hw],'xticklabel',[0 300 600],'yticklabel',[0 300 600]);
    hold on;
    i = 0;
    j = 0;
    for t = tp{id_out-1} % 1:length(t_mid)
        h2 = plot(100,0);
        h3 = plot(100,0);
        ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
        h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'b.');
        if ~isnan(x_centre(t)) && max(abs(R.grid.quick.centre(:,t))) < 31.5
            x_tmp = x_centre(t);
            y_tmp = y_centre(t);
            %             [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
            %             spikingn = unique(spikingn);
            %             all = find(lattice_nD_find_dist(Lattice,hw,round(x_tmp),round(y_tmp)) <= width(t));
            %             incircle = sum(ismember(spikingn,all));
            if 1 % incircle/length(all) >= 0.07 || incircle > 23
                if i == 0
                    x0 = x_tmp;
                    y0 = y_tmp;
                end
                i = i + 1;
                r_cos = x_tmp+width(t)*cos(ang);
                r_sin = y_tmp+width(t)*sin(ang);
                h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'k', ...
                    r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'k',...
                    r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'k',...
                    r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'k',...
                    r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'k', ...
                    r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'k', ...
                    r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'k', ...
                    r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'k', ...
                    r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'k','LineWidth',1.5);
                h2 = plot( x_tmp, y_tmp, 'r>'); % , 'MarkerSize', 8);
                if sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18 % Distance_xy(x_tmp,y_tmp,x0,y0,fw) < 18
                    plot([x0,x_tmp],[y0,y_tmp],'r','LineWidth',1.5)
                    j = 0;
                else
                    j = 1;
                end
                x0 = x_tmp;
                y0 = y_tmp;
            end
        end        
        if t ~= tp{id_out-1}(end)
            delete(h1);
            delete(h2);
            delete(h3);
        else
            text(-0.28,1,str2{id_out-1},'Units', 'Normalized','FontSize',12)
        end
        if j == 1
            delete(findobj(gca,'Type','line','Color','r'));
        end
%         if id_out > 2
%             ts = str{id_out-2};
%             title(ts);
%         end
    end
end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc GammaPaperFig4V3 % this is the trick!!