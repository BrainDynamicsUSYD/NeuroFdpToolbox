function d = Velocity(R)
% calculate modified spike pattern move velocity
% jump distacne is velocity per ms
coe = 2; % 2
M = 20; % 20
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
t_mid = R.grid.t_mid;
centre = R.grid.quick.centre;
width = R.grid.quick.radius;
for t = 1:length(t_mid)
    if ~isnan(centre(1,t))
        x_tmp = centre(1,t);
        y_tmp = centre(2,t);
        
        % adding modification algorithm as criteria for proper spike pattern
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        spikingn = unique(spikingn);
        all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        incircle = sum(ismember(spikingn,all));
        if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
        else
            centre(:,t) = NaN;
        end
    end
end
d = Distance_xy(centre(1,1:end-1),centre(2,1:end-1),centre(1,2:end),centre(2,2:end),2*hw+1);
end