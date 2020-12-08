function VisualizeGridFiringCentre4(R,mode)
% adapt from function VisualizeGridFiringCentre.m
% use complex estimator to decide pattern
% + raster plot
% + wavelet spectrogram
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
fig = figure('Name','Vis','color','w','NumberTitle','off');
ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

% Butterworth filter
order = 4; % 4th order
lowFreq = 30; % broad band (1-1000 Hz)
hiFreq = 80;
dt = R.dt;
fs = 1/(1e-3); % *dt); % sampling frequency (Hz)
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

% LFP = filter(b,a,R.LFP.LFP{1}(5,:));
% fc = centfrq('cmor1.5-1');
% scalerange = fc./([lowFreq hiFreq]*(1/fs));
% scales = scalerange(end):0.5:scalerange(1);
% pseudoFreq = scal2frq(scales,'cmor1.5-1',1/fs);

% FR = R.Analysis.Hz_t{1};
% FR = R.num_spikes{1}/3969*1e4;

FR = R.grid.num_spikes_win/3969*200; % Hz 
FRcertain = filter(b,a,FR);

i = 0;
j = 0;
mark = 1;
coe = 2;
M = 20;

% R = load('0001-201806261601-18754_in_1529993007138_out_RYG.mat', 'LFP');
% tall = length(ripple_power_grid);
% mm = minmax(ripple_power_grid(:)');
% Y = prctile(ripple_power_grid(:),96);
% mm(1) = Y;

% [R] = GetBurst2(R);
% Pattern2D = reshape((R.LFP.GammaBurstEvent.is_burst).*R.LFP.LFP_gamma_hilbert_abs,20,20,[]);
% mm = minmax(Pattern2D(:)');

% vidObj = VideoWriter('0005Dynamics.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 10; % number of frames to display per second
% open(vidObj);
% set(gcf,'Position',[680 284 917 677])
%%
mark = 1;
for t = 1:length(t_mid) % 1:length(t_mid)
%     subplot(1,2,1)
    subplot(3,2,[1 3 5])
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
        
        % adding modification algorithm as criteria for proper spike pattern
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        spikingn = unique(spikingn);
        all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        incircle = sum(ismember(spikingn,all));
        if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
            
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
%             if mark == 1
%                 plot([-31,-31;31,31],[-21,-11;-21,-11],'k')
%                 mark = 0;
%             end
            
            if Distance_xy(x0,y0,x_tmp,y_tmp,63) < 1 % sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18
                plot([x0,x_tmp],[y0,y_tmp],'g')
                j = 0;
            else
                j = 1;
            end
            
        end
    end
    x0 = x_tmp;
    y0 = y_tmp;
    %     drawnow
%     F = getframe(fig);
%     writeVideo(vidObj,F.cdata);

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
    
%     subplot(1,2,2)
%     imagesc(flipud(ripple_power_grid(:,:,t_mid(t)-2e4+1)'),mm)
%     imagesc(flipud(Pattern2D(:,:,t_mid(t))'),mm)
%     imagesc(Pattern2D(:,:,t_mid(t)),mm)
            
%     subplot(3,2,2) % raster plot
%     hold on;
%     if t == 1 || mod(time_ms,1e3) < 1 % 5e2
%         seg = 1e4*floor(time_ms/1e3)+1:1e4*(floor(time_ms/1e3)+1); % 5e3
%         if max(seg) > R.step_tot
%             a = 0;
%         else
%             a = 1;
%         end
%         RasterPlotYL(R,1,1,[],ceil(dt*seg(1)):ceil(dt*seg(end))+a)
%     end

    subplot(3,2,2) % raster plot
    hold on;
    if t == 1 || mod(time_ms,1e3) < 1 % 5e2
        seg = 1e3*floor(time_ms/1e3)+1:1e3*(floor(time_ms/1e3)+1); % 5e3
        [wt,f,coi] = cwt(FRcertain(seg),'morse',fs,'TimeBandwidth',60);
        ind1 = find(f >= lowFreq & f <= hiFreq);
        ind2 = find(coi<30);
        f = f(ind1);
        wt = wt(ind1,:);
        temporal = seg; % *dt;
        temporal = temporal(ind2);
        fr = FRcertain(seg);
        fr = fr(ind2);
        if max(seg) > R.step_tot
            a = 0;
        else
            a = 1;
        end
        RasterPlotYL(R,1,1,[],ceil(temporal(1)):ceil(temporal(end))+a)
    end
    
    subplot(3,2,4) % temporal LFP/rate signal
    hold on;

%     if t == 1 || mod(time_ms,5e2) < 1
%         hold off
%         plot(dt*seg,R.LFP.LFP{1}(5,seg),'color',[0.8 0.8 0.8])
%         hold on;
%         plot(dt*seg,LFP(seg),'k')
%         legend('raw LFP','SWR')
%         %         text(-0.1,1.12,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%         ylabel('LFP(a.u.)')
%     end

    if t == 1 || mod(time_ms,1e3) < 1
        hold off
        plot(temporal,abs(fr))
        xlim(temporal([1 end]));
        ylabel('Population rate(Hz)')
    end
    
    subplot(3,2,6);
    hold on;
    if t == 1 || mod(time_ms,1e3) < 1 % 5e2
        hold off
        CData = abs(wt(end:-1:1,ind2));
        Y = prctile(CData(:),80);
        CData(CData < Y) = 0;
        uimagesc(temporal,f(end:-1:1),CData)
        
%         uimagesc(temporal,f(end:-1:1),abs(wt(end:-1:1,ind2)))
%         coeffs_tmp = abs(cwt(LFP(seg),scales,'cmor1.5-1'))'; % this had narrow width of scales
%         CData = transpose(coeffs_tmp); % coeffs_tmp'
%         uimagesc(dt*seg,pseudoFreq(end:-1:1),CData(end:-1:1,:));

        colorbar('off');
        
        %         xlim([floor(time_ms/5e2) floor(time_ms/5e2)+1]*5e2) % unit:ms
        %     xlim([1*(seg-1) 1*seg]);
        
        ylim([lowFreq hiFreq]);
        ylabel('Frequency(Hz)');
        set(gca,'YDir','normal')
    end
    hold on;
    delete(findobj(gca,'Type','line','Color','k'));
    xlabel('Time(ms)')
    
%     xlim([floor(time_ms/5e2) floor(time_ms/5e2)+1]*5e2) % unit:ms

    plot([time_ms time_ms],[0 1e5],'k-');
    
    %     ylim([0 ymax])
end
% close(gcf);
% close(vidObj);
end
