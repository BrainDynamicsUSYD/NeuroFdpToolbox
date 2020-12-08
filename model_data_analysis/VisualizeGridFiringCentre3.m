function VisualizeGridFiringCentre3(R,mode,IndC,LFP_centre_x,LFP_centre_y)
% adapt from function VisualizeGridFiringCentre.m
% use complex estimator to decide pattern
% + theta/beta/gamma time-vary power curve
close all;
clc;
% mode = 'bayesian';

hw = 31;
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


ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

% Butterworth filter
order = 4; % 4th order
lowFreq = [4 20 30]; % theta/beta/gamma
hiFreq = [7 30 80];
fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
LFP = R.LFP.LFP{1};
ca = cell(1,length(IndC)*3);
LFP_power = zeros(3*length(IndC),length(LFP)); % [t;t;b;b;g;g]
for k = 1:3
    Wn = [lowFreq(k) hiFreq(k)]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
    for j = 1:length(IndC)
        [~,indE] = sort(Distance_xy(LFP_centre_x,LFP_centre_y,Lattice(IndC(j),1),Lattice(IndC(j),2),2*hw));
        LFP_tbg = filter(b,a,LFP(indE(1),:));
        LFP_power((k-1)*length(IndC)+j,:) = abs(hilbert(LFP_tbg)).^2/(hiFreq(k)-lowFreq(k));
        ca{(k-1)*length(IndC)+j} = sprintf('Electrode No.%d',indE(1));
    end
end
colorv = ['g' 'b' 'r'];
i = 0;
post_dist = lattice_nD_find_dist(Lattice,hw,513); 
[r,~] = sort(post_dist);
r = r(R.ExplVar.local_population);
clear post_dist
mark = 1;
for t = 1:length(t_mid)
    subplot(2,1,1)
    axis equal;
    box on;
    set(gca,'xtick',[],'ytick',[]);
    
    xlim([-hw hw]);
    ylim([-hw hw]);
    hold on;
    h2 = plot(100,0);
    h3 = plot(100,0);
    
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    if ~isnan(x_centre(t))
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        
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
        if nargin >= 3 && mark == 1
            %             hold on;
            %             h4 = viscircles(Lattice(hotpot,:),5.66,'Color','y');
            for k = 1:length(IndC)
                r_cos = Lattice(IndC(k),1) + r*cos(ang); % N400:11.2, N100:5.66 N200:8.1
                r_sin = Lattice(IndC(k),2) + r*sin(ang);
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
        
        if Distance_xy(x0,y0,x_tmp,y_tmp,63) < 10 % sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18
            plot([x0,x_tmp],[y0,y_tmp],'g')
            j = 0;
        else
            j = 1;
        end
        
    end
    x0 = x_tmp;
    y0 = y_tmp;
    %     drawnow
    pause(0.05);
    delete(h1);
    delete(h2);
    delete(h3);
    
    if j == 1
        delete(findobj(gca,'Type','line','Color','g'));
    end
    time_ms = t_mid(t)*0.1;
    ts = sprintf('time = %8.1f ms',time_ms);
    title(ts);
    
    subplot(2,1,2);
    hold on;
    if t == 1 || mod(time_ms,5e2) < 1
        seg = 5e3*floor(time_ms/5e2)+1:5e3*(floor(time_ms/5e2)+1);
        for k = 1:3
            for m = 1:length(IndC)
                plot(seg*0.1,LFP_power((k-1)*length(IndC)+m,seg),colorv(k),'LineWidth',0.75*m); % theta green;beta blue;gamma red
%                 hold on;
            end
        end
        ymax = max(max(LFP_power(:,seg)));        
        legend(ca);
    end
    delete(findobj(gca,'Type','line','Color','k'));
    xlabel('Time(ms)')
    ylabel('Power(mV^2/Hz)')
    xlim([floor(time_ms/5e2) floor(time_ms/5e2)+1]*5e2) % unit:ms
    plot([time_ms time_ms],[0 1e5],'k-');
    ylim([0 ymax])
    
    
    %     LFP_power_tbg(:,1) = LFP_power_tbg(:,2);
    %     for k = 1:3
    %         for m = 1:length(IndC)
    %             LFP_power_tbg((k-1)*length(IndC)+m,2) = LFP_power((k-1)*length(IndC)+m,t_mid(t));
    %             plot([t_mid(t-1) t_mid(t)]*0.1,LFP_power_tbg((k-1)*length(IndC)+m,:),colorv(k),'LineWidth',0.75*m); % theta green;beta blue;gamma red
    % %             hold on;
    %         end
    %     end
    %     legend(ca); % ([1 3 5]));
    %     xlim([floor(time_ms/2e2) floor(time_ms/2e2)+1]*2e2) % unit:ms
end

end
