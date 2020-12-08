function TuningParameterForPattern(varargin)
%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

%%% basic setting %%%
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 0; % number of trials - 1
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 1; % number of calculating trails
Opt.testMode = 1;
[Lattice2, ~] = lattice_nD(2, 1.5); % creat 4*4
bin = 300; % 30ms
nonoverlap = bin; % unit: 0.1ms
start = 10006; % R.grid.t_mid = 26:10:99966
No = 1; % No. of Electrode
dc = 0;

%%% Loop number for PBS array job %%%
loop_num = 0;

ts = [];
num = [];
hw = 31;
[Lattice, ~] = lattice_nD(2, hw);
turingt = ((ones(3,1)*13*[0:9])' + [6:8])';
turingt = turingt(:);
for r = turingt'
    for i = [3]
        for coe = [1.0]
            for th1 = [0.05 0.07 0.1]
                for th2 = [20 23 25]
                    
                    % For PBS array job
                    loop_num = loop_num + 1;
                    if nargin ~= 0
                        PBS_ARRAYID = varargin{1};
                        if loop_num ~=  PBS_ARRAYID
                            continue;
                        end
                    end
                    
                    fprintf('Loading RYG.mat file %s...', files{r});
                    R = load(files{r});
                    R = get_grid_firing_centre(R);
                    t_mid = R.grid.t_mid;
                    x_centre = R.grid.quick.centre(1,:);
                    y_centre = R.grid.quick.centre(2,:);
                    width = R.grid.quick.radius;
                    [no,steps] = size(R.LFP.LFP_broad);
                    LFP_certain = zeros(no,steps);
                    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
                    
                    for t = 1:length(t_mid)
                        if ~isnan(x_centre(t)) && max(abs(R.grid.quick.centre(:,t))) < 31.5
                            x_tmp = x_centre(t);
                            y_tmp = y_centre(t);
                            [spikingn,~] = find(R.spike_hist{1}(:,(t_mid(t)-25):(t_mid(t)+24))); % because win_len = 50 in get_grid_firing_centre.m
                            spikingn = unique(spikingn);
                            all = find(lattice_nD_find_dist(Lattice,hw,x_tmp,y_tmp) <= width(t)); % 81~1257
                            incircle = sum(ismember(spikingn,all));
                            %             if incircle/length(all) >= 0.07 && incircle > 23 % old setting: >=0.07 || >23
                            switch i
                                case 1
                                    if incircle/length(spikingn) > coe*pi*width(t)^2/(2*hw+1)^2 && incircle > th
                                        ts = [ts t_mid(t)*1e-4]; % s
                                        num = [num incircle];
                                    end
                                case 2
                                    if incircle/length(spikingn) > coe*pi*width(t)^2/(2*hw+1)^2 || incircle > th
                                        ts = [ts t_mid(t)*1e-4]; % s
                                        num = [num incircle];
                                    end
                                case 3
                                    if incircle/length(spikingn) > coe*pi*width(t)^2/(2*hw+1)^2 && ( incircle/length(all) > th1 || incircle > th2 )
                                        ts = [ts t_mid(t)*1e-4]; % s
                                        num = [num incircle];
                                    end
                            end
                        end
                    end
                    pt = 1e4*ts; % P.ts: second  pt: time step
                    
                    % set seed
                    date_now = datestr(now,'yyyymmddHHMM-SSFFF');
                    scan_temp = textscan(date_now,'%s','Delimiter','-');
                    rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
                    alg='twister';
                    rng(rand_seed,alg); %  Note that this effect is global!
                    
                    % Butterworth filter
                    order = 4; % 4th order
                    lowFreq_br = 5;
                    hiFreq_br = 100;
                    Wn = [lowFreq_br hiFreq_br]/(fs/2);
                    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
                    for j = 1:no
                        LFP_certain(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
                    end
                    Hilbert = zeros(size(LFP_certain));
                    for j = 1:no
                        Hilbert(j,:) = hilbert(LFP_certain(j,:));
                    end
                    send = No;
                    d = lattice_nD_find_dist(Lattice2,1.5,send);
                    receive = find(d>0)'; % d > 1.5
                    for tt = start:nonoverlap:steps-2*bin
                        period = tt:(tt+bin-1);
                        ind = find(ismember(period,pt));
                        if ~isempty(ind)
                            t = period(ind(1));
                        else
                            continue
                        end
                        for k = 1:length(receive)
                            nte = zeros(1,bin);
                            ntesh = zeros(1,bin);
                            X = repmat(angle(Hilbert(send,t:(t+bin-1)))+pi,4,1);
                            Y = repmat(angle(Hilbert(receive(k),t:(t+bin-1)))+pi,4,1);
                            x = X(1,:);
                            y = Y(1,:);
                            Xsh = repmat(x(randperm(bin)),4,1);
                            Ysh = repmat(y(randperm(bin)),4,1);
                            for j = -(1:bin)
                                Opt.tauy = j;
                                try
                                    [NTE] = transferentropy(X',Y',Opt,'NTE');
                                    NTE = NTE(NTE ~= Inf);
                                    NTE = NTE(NTE ~= -Inf);
                                    [NTEsh] = transferentropy(Xsh',Ysh',Opt,'TE','NTE');
                                    NTEsh = NTEsh(NTEsh ~= Inf);
                                    NTEsh = NTEsh(NTEsh ~= -Inf);
                                    ntesh(abs(j)) = nanmean(NTEsh);
                                    if nanmean(NTE) > 0
                                        nte(abs(j)) = nanmean(NTE);
                                    end
                                catch
                                end
                            end
                            if max(nte) > max([max(ntesh),0])
                                dc = dc + 1;
                            end
                        end
                    end
                    save(['T',sprintf('%04g',r),'ratio',sprintf('%0.2f',th1),'th',sprintf('%d',th2),'NTE.mat'],'dc')
                end
            end
        end
    end
end
end