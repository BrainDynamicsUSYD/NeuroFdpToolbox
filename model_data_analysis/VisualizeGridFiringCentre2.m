function VisualizeGridFiringCentre2(R,mode,hotpot,duration)
% adapt from function visualize_grid_and_SWR.m
% use complex estimator to decide pattern
% optional: adding local hotpot
% close all;
clc;
% mode = 'bayesian';

hw = 31;
fw = 2*hw+1;
t_mid = R.grid.t_mid;
% ind_a_vec = R.grid.ind_ab(1,:);
% ind_b_vec = R.grid.ind_ab(2,:);

switch mode
    case 'bayesian'
        x_centre = R.grid.bayes.centre(1,:);
        y_centre = R.grid.bayes.centre(2,:);
        width = R.grid.bayes.radius;
    case 'bayesian2'
        x_centre = R.grid.bayes.centre(1,:);
        y_centre = R.grid.bayes.centre(2,:);
        width = R.grid.bayes.radius;
        for i = 1:length(t_mid)
            if R.grid.bayes.bayes_factor_ln(i) < log(100)
                x_centre(i) = NaN;
                y_centre(i) = NaN;
                width(i) = NaN;
            end
        end
    case 'quick'
        x_centre = R.grid.quick.centre(1,:);
        y_centre = R.grid.quick.centre(2,:);
        width = R.grid.quick.radius;
end

[Lattice, ~] = lattice_nD(2, hw);
% x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
% y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);

% figure('Name','Vis','color','w','NumberTitle','off');
axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);

xlim([-hw hw]);
ylim([-hw hw]);
hold on;

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

% i = 0;
mark = 1;
post_dist = lattice_nD_find_dist(Lattice,hw,513);
[r,~] = sort(post_dist);
r = r(R.ExplVar.local_population);
clear post_dist

% coe = 2;
% M = 20;
% duration = 2.4e3:5e3; % 1:1.9e3;
for t = duration % 2.4e3:5e3 % 2e3:2.35e3 % 1:1.9e3 % 2.4e3:length(t_mid) % 1:length(t_mid)
    h2 = plot(100,0);
    %     h3 = plot(100,0);
    
    %     ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    %     h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    if ~isnan(x_centre(t))
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        
        % adding modification algorithm as criteria for proper spike pattern
%         [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
%         spikingn = unique(spikingn);
%         all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
%         incircle = sum(ismember(spikingn,all));
%         if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
            
            %         if i == 0
            %             x0 = x_tmp;
            %             y0 = y_tmp;
            %         end
            %         i = i + 1;
            
            %         r_cos = x_tmp+width(t)*cos(ang);
            %         r_sin = y_tmp+width(t)*sin(ang);
            %         h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
            %             r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
            %             r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
            %             r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
            %             r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
            %             r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
            %             r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
            %             r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
            %             r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
            
            h2 = plot( x_tmp, y_tmp, 'r.', 'MarkerSize', 8);
            if mark == 1
                for k = 1:length(hotpot)
                    r_cos = Lattice(hotpot(k),1) + r*cos(ang); % N400:11.2, N100:5.66 N200:8.1
                    r_sin = Lattice(hotpot(k),2) + r*sin(ang);
                    plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'k', ...
                        r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'k',...
                        r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'k',...
                        r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'k',...
                        r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'k', ...
                        r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'k', ...
                        r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'k', ...
                        r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'k', ...
                        r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'k');
                end
                mark = 0;
            end
%         end
        
        %         if sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18
        %             plot([x0,x_tmp],[y0,y_tmp],'g')
        %             j = 0;
        %         else
        %             j = 1;
        %         end
        
    end
    %     x0 = x_tmp;
    %     y0 = y_tmp;
    %     pause(0.05);
    %     delete(h1);
    %     delete(h2);
    %     delete(h3);
    
    %     if j == 1
    %         delete(findobj(gca,'Type','line','Color','g'));
    %     end
    %     if t_mid(t)*0.1 > 50e3
    %     ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);
    %     title(ts);
    %     break
    %     end
end
ts = sprintf('time = %8.1f ms', (t_mid(duration(end))-t_mid(duration(1)))*0.1);
title(ts);
end
