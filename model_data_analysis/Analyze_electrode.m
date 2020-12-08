function Analyze_electrode(R)
dt = R.dt;
t_mid = R.grid.raw.t_mid;
Electodes = zeros(2785,6);
vec_time = [];
local = [1 2 4 5 6 8 13 14 16;
    1 2 3 5 6 7 13 14 15;
    2 3 4 6 7 8 14 15 16;
    1 3 4 5 7 8 13 15 16;
    1 2 4 5 6 8 9 10 12;
    1 2 3 5 6 7 9 10 11;
    2 3 4 6 7 8 10 11 12;
    1 3 4 5 7 8 9 11 12;
    5 6 8 9 10 12 13 14 16;
    5 6 7 9 10 11 13 14 15;
    6 7 8 10 11 12 14 15 16;
    5 7 8 9 11 12 13 15 16;
    1 2 4 9 10 12 13 14 16;
    1 2 3 9 10 11 13 14 15;
    2 3 4 10 11 12 14 15 16;
    1 3 4 9 11 12 13 15 16];
i = 1;
m = 1;
text1 = 'Start time:%0.5g ms  End time:%0.5g ms\n';

for tt = unique(sort([1:10:t_mid(end) t_mid(:)']))
    
    % ripple detection
    detceded = sum(R.LFP.ripple_event.is_SWR(:, tt), 2) > 0;
    ind = find(detceded);
    if ~isempty(ind)
        Electodes(i,1:length(ind)) = ind;
        i = i + 1;
        vec_time = [vec_time tt*dt];
    end
end
while m < 2785
    around = [1:16];
    ts = [];
    while sum(ismember(Electodes(m,:),around)) ~= 0
        ts = [ts vec_time(m)];
        ind_elec = Electodes(m,:);
        ind_elec = ind_elec(ind_elec ~= 0);
        around = unique(local(ind_elec,:));
        m = m + 1;
    end
    if length(ts) > 1
        fprintf(text1,ts(1),ts(end));
    end
end
end