function [Jum,Ste,Rad,Dur,Spi] = CollectDurationJump(R,mode,thresh)
% jump size; steady size; radius; duration; step_spikes
hw = 31;
t_mid = R.grid.t_mid;
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
i = 0;
Jum = [];
Ste = [];
Rad = [];
Spi = cell(1,2);
tind = zeros(1,length(t_mid));
coe = 2;
M = 20;
for t = 1:length(t_mid)
    if ~isnan(x_centre(t))
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        % adding modification algorithm as criteria for proper spike pattern
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        spikingn = unique(spikingn);
        all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        incircle = sum(ismember(spikingn,all));
        if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
            if y_tmp > -11 || y_tmp < -21
                i = 0;
                continue
            end
            if i == 0
                x0 = x_tmp;
                y0 = y_tmp;
                i = i + 1;
                continue
            end
            d = Distance_xy(x0,y0,x_tmp,y_tmp,63);
            if d < thresh % sqrt((x_tmp - x0)^2 + (y_tmp - y0)^2) < 18
                Ste = [Ste d];
                tind(t) = 1;
            else
                Jum = [Jum d];
            end
            Rad = [Rad width(t)];            
            Spi{1} = [Spi{1} d];
            Spi{2} = [Spi{2} sum(sum(full(R.spike_hist{1}(630:1260,t-35:t+24))))];
        end
    end
    x0 = x_tmp;
    y0 = y_tmp;
end
[~,Dur,~,~,~] = seq_postprocess(tind,1);
end