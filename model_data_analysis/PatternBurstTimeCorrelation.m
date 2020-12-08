function [a,b,c,pr] = PatternBurstTimeCorrelation(R)
hw = 31;
t_mid = R.grid.t_mid;

x_centre = R.grid.quick.centre(1,:);
y_centre = R.grid.quick.centre(2,:);
width = R.grid.quick.radius;

[Lattice, ~] = lattice_nD(2, hw);


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

FR = R.grid.num_spikes_win/3969*200; % Hz
FRcertain = filter(b,a,FR);

coe = 2;
M = 20;
time_ms1 = [];
for t = 1:length(t_mid)-1e3 % 1:length(t_mid)
    if ~isnan(x_centre(t))
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        % adding modification algorithm as criteria for proper spike pattern
        [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24)));
        spikingn = unique(spikingn);
        all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t));
        incircle = sum(ismember(spikingn,all));
        if (incircle/(pi*width(t)^2) >= coe*length(spikingn)/((2*hw)^2)) && (incircle >= M)
            time_ms1 = [time_ms1 t_mid(t)*0.1];
        end
    end
end
time_ms2 = [];
for i = 1:1e3:length(FR)-999
    seg = 1e3*floor(i/1e3)+1:1e3*(floor(i/1e3)+1); % 5e3
    [wt,f,coi] = cwt(FRcertain(seg),'morse',fs,'TimeBandwidth',4);
    ind1 = find(f >= lowFreq & f <= hiFreq);
    ind2 = find(coi<30);
    period = seg(1)-1+ind2; % ms
    f = f(ind1);
    wt = wt(ind1,:);
    temporal = seg; % *dt;
    temporal = temporal(ind2);
    fr = FRcertain(seg);
    fr = fr(ind2);
    CData = abs(wt(end:-1:1,ind2));
    Y = prctile(CData(:),80);
    CData(CData < Y) = 0;
    CData = sum(CData);
    time_ms2 = [time_ms2 period(CData > 0)']; % ms
end
% r = xcorr(time_ms1,time_ms2);
% plot(time_ms1,time_ms1,'ro',time_ms2,time_ms2,'k*')
x = ismember(1:(length(FR)-1e3),round(time_ms1));
y = ismember(1:(length(FR)-1e3),time_ms2);
% [r,lags] = xcorr(x,y);
% plot(lags,r)

a = length(find(x==y));
b = length(find(x));
c = length(find(y));
pr = a/(length(FR)-1e3);
end