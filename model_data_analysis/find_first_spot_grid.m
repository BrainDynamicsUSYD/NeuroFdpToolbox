function [Ind,tid] = find_first_spot_grid(R)
ind = find(~isnan(R.grid.quick.centre(1,:)));
Ind = [];
for i = 2:length(ind)
    if ind(i) > ind(i-1) + 1
        Ind = [Ind ind(i)];
    end
end
tp = 0.1*R.grid.t_mid(Ind); % ms (time point)
tid = tp(2:end) - tp(1:end-1); % (time interval distribution)
% histogram(tid,20)
% xlabel('ms')
end