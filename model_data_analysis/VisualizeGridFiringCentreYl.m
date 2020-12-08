function VisualizeGridFiringCentreYl(R,mode)
% adapt from function visualize_grid_and_SWR.m
% use complex restimator + my modification to decide pattern
close all;
clc;
% mode = 'bayesian';
try
    hw = R.ExplVar.hw;
catch
    hw = 31;
end
fw = 2*hw+1;
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

[Lattice, ~] = lattice_nD(2, hw);
x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);

figure('Name','Vis','color','w','NumberTitle','off');
axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);

xlim([-hw hw]);
ylim([-hw hw]);
hold on;

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

i = 0;
j = 0;
for t = 750:810 % 1:length(t_mid)
    h2 = plot(100,0);
    h3 = plot(100,0);
    
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    if ~isnan(x_centre(t)) && max(abs(R.grid.quick.centre(:,t))) < 31.5
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        spikingn = unique(spikingn);
        all = find(lattice_nD_find_dist(Lattice,hw,round(x_tmp),round(y_tmp)) <= width(t));
        incircle = sum(ismember(spikingn,all));
        if incircle/length(all) >= 0.07 || incircle > 23
            if i == 0
                x0 = x_tmp;
                y0 = y_tmp;
            end
            i = i + 1;
            
            r_cos = x_tmp+width(t)*cos(ang);
            r_sin = y_tmp+width(t)*sin(ang);
            h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
                r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
                r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
                r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
                r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
                r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
                r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
                r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
                r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
            
            h2 = plot( x_tmp, y_tmp, 'r>', 'MarkerSize', 8);
            
            if sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18 % Distance_xy(x_tmp,y_tmp,x0,y0,fw) < 18
                plot([x0,x_tmp],[y0,y_tmp],'g')
                j = 0;
            else
                j = 1;
            end
            x0 = x_tmp;
            y0 = y_tmp;
        end
    end
    text(-0.07,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%     if t_mid(t) > 9000
%         break
%     end
    pause(0.05);
    delete(h1);
% if t_mid(t) < 7500
    delete(h2);
% end
    delete(h3);
    
    if j == 1
        delete(findobj(gca,'Type','line','Color','g'));
    end
    
    ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);  
%     ts = sprintf('time from 750 ms to 900 ms'); 
    title(ts);    
end

end
