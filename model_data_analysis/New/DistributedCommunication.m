function DistributedCommunication(varargin)
% method "1X" means working on spike trains
%        11 : random select spots as send and receive;
%             tauy needed in calculating transfer entropy is dependent on
%             distance between two sites;
%             after start,calculate te in bin every bin-nonoverlap among
%             spots-1 pairs;
%             balance information amount and distributed ability to get dc
%        12 : send is the location of pattern;
%             others are nearly all the same
%        "2X" means working LFP phase
%        21 : select electrode No as send and others as receive;
%        22 : send is the electrode which has pattern;
%             others are all the same
% Both methods work on bins which MUST contain patterns.
% I haven't checked menthod 11 & 12.
tic;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

dir_strut2 = dir('Pattern*.mat');
files2 = cell(1,num_files);
for id_out = 1:num_files
    files2{id_out} = dir_strut2(id_out).name;
end

%%% basic setting %%%
hw = 31;
spots = 7;
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 0; % number of trials - 1
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 1; % number of calculating trails
Opt.testMode = 1;
dmin = 0;
DC = [];
METHOD = [];
[Lattice1, ~] = lattice_nD(2, hw);
[Lattice2, ~] = lattice_nD(2, 1.5); % creat 4*4
LFP_centre_x = linspace(-hw, hw, 9);
LFP_centre_y = linspace(-hw, hw, 9);
LFP_centre_x = LFP_centre_x(2:2:8);
LFP_centre_y = LFP_centre_y(2:2:8);
[LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
LFP_centre_x = LFP_centre_x(:);
LFP_centre_y = LFP_centre_y(:);

%%% Loop number for PBS array job %%%
loop_num = 0;

%%%%%%%%%%%%%%% variables %%%%%%%%%%%%%%%%%%
% v = [0.64 0.61 0.60 0.61 0.71 1.20 1.63 1.96 2.05 2.05 1.97 1.85 1.88]; % grid_d/ms:v*10um/ms  um/ms = mm/s

bin = 300; % 30ms
nonoverlap = bin; % unit: 0.1ms
start = 10006; % R.grid.t_mid = 26:10:99966

%%% mode 1:spikes %%%

%%% mode 2:LFP %%%
No = 1; % No. of Electrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:num_files
     
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    P = load(files2{i});
    R = get_grid_firing_centre(R);
    [no, steps] = size(R.LFP.LFP_broad);
    pt = 10*P.ts;
    px = P.x;
    py = P.y;
    pr = P.r;
    v = nanmean(P.v); % d_grid/ms
    M = [21];
    for method = M % [11 12]
        dc = 0;
        
        % set seed
        date_now = datestr(now,'yyyymmddHHMM-SSFFF');
        scan_temp = textscan(date_now,'%s','Delimiter','-');
        rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
        alg='twister';
        rng(rand_seed,alg); %  Note that this effect is global!
        
        switch method
            case 11
                ntesh = zeros(1,spots-1);
                while dmin < 20
                    b = randperm(length(Lattice1),spots-1);
                    for j = 1:spots-1
                        dist = lattice_nD_find_dist(Lattice1(b,:),hw,Lattice1(b(j),1),Lattice1(b(j),2));
                        dmin = min(dist(dist>0));
                        if dmin < 20
                            break
                        end
                    end
                end
                send = randperm(length(Lattice1),1);
                receive = b;
                tauy = -round(Distance_xy(Lattice1(receive,1),Lattice1(receive,2),Lattice1(send,1),Lattice1(send,2),hw*2 + 1)/v*10); % 0.1ms
                for tt = start:nonoverlap:steps-bin
                    period = tt:(tt+bin-1);
                    ind = find(ismember(period,pt));
                    if ~isempty(ind)
                        t = period(ind(1));
                    else
                        continue
                    end
                    for r = 1:spots-1
                        Opt.tauy = tauy(r); % 0.1ms
                        try
                            [NTEsh] = transferentropy(R.spike_hist{1}(FindNeurons(Lattice1,hw,10,send)',t:(t+bin-1))',R.spike_hist{1}(FindNeurons(Lattice1,hw,10,receive(r))',t:(t+bin-1))',Opt,'NTEsh');
                            NTEsh = NTEsh(NTEsh ~= Inf);
                            NTEsh = NTEsh(NTEsh ~= -Inf);
                            if nanmean(NTEsh) > 0
                                ntesh(r) = nanmean(NTEsh);
                            end
                        catch
                        end
                    end
                    dte = sum(ntesh);
                    %                     numte = length(find(ntesh));
                    %                     switch numte
                    %                         case {0,1}
                    %                             dte = 0;
                    %                         otherwise
                    %                             dte = numte*sum(ntesh);
                    %                     end
                    dc = dc + dte;
                end
                METHOD = [METHOD method];
                DC = [DC dc];
            case 12
                ntesh = zeros(1,spots-1);
                while dmin < 20
                    b = randperm(length(Lattice1),spots-1);
                    for j = 1:spots-1
                        dist = lattice_nD_find_dist(Lattice1(b,:),hw,Lattice1(b(j),1),Lattice1(b(j),2));
                        dmin = min(dist(dist>0));
                        if dmin < 20
                            break
                        end
                    end
                end
                receive = b;
                for tt = start:nonoverlap:steps-bin
                    period = tt:(tt+bin-1);
                    ind = find(ismember(period,pt));
                    if ~isempty(ind)
                        t = period(ind(1));
                    else
                        continue
                    end
                    cs(1) = round(px(pt == t));
                    cs(2) = round(py(pt == t));
                    if sum(isnan(cs)) > 0 || max(abs(cs)) > hw
                        continue
                    else
                        send = IndexLattice(Lattice1,cs(1),cs(2));
                    end
                    tauy = -round(Distance_xy(Lattice1(receive,1),Lattice1(receive,2),Lattice1(send,1),Lattice1(send,2),hw*2 + 1)/v*10); % 0.1ms
                    [spikingn,~] = find(R.spike_hist{1}(:,(t-25):(t+24)));
                    spikingn = unique(spikingn);
                    all = find(lattice_nD_find_dist(Lattice1,hw,px(pt == t),py(pt == t)) <= pr(pt == t));
                    incircle = spikingn(ismember(spikingn,all));
                    for r = 1:spots-1
                        Opt.tauy = tauy(r); % 0.1ms
                        try
                            [NTEsh] = transferentropy(R.spike_hist{1}(incircle(randperm(length(incircle),10)),t:(t+bin-1))',R.spike_hist{1}(FindNeurons(Lattice1,hw,10,receive(r))',t:(t+bin-1))',Opt,'NTEsh');
                            NTEsh = NTEsh(NTEsh ~= Inf);
                            NTEsh = NTEsh(NTEsh ~= -Inf);
                            if nanmean(NTEsh) > 0
                                ntesh(r) = nanmean(NTEsh);
                            end
                        catch
                        end
                    end
                    dte = sum(ntesh);
                    %                     numte = length(find(ntesh));
                    %                     switch numte
                    %                         case {0,1}
                    %                             dte = 0;
                    %                         otherwise
                    %                             dte = numte*sum(ntesh);
                    %                     end
                    dc = dc + dte;
                end
                METHOD = [METHOD method];
                DC = [DC dc];
            case 21
                try
                    LFP_certain = R.LFP.LFP_broad;
                catch
                    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
                    
                    % Butterworth filter
                    order = 4; % 4th order
                    lowFreq_br = 5;
                    hiFreq_br = 100;
                    Wn = [lowFreq_br hiFreq_br]/(fs/2);
                    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
                    for j = 1:no
                        LFP_certain(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
                    end
                end
                Hilbert = zeros(size(LFP_certain));
                for j = 1:no
                    Hilbert(j,:) = hilbert(LFP_certain(j,:));
                end
                send = No;
                d = lattice_nD_find_dist(Lattice2,1.5,send);
                receive = find(d>0)';
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
                METHOD = [METHOD method];
                DC = [DC dc];
            case 22
                try
                    LFP_certain = R.LFP.LFP_gamma;
                catch
                    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
                    
                    % Butterworth filter
                    order = 4; % 4th order
                    lowFreq_br = 5;
                    hiFreq_br = 100;
                    Wn = [lowFreq_br hiFreq_br]/(fs/2);
                    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
                    for j = 1:no
                        LFP_certain(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
                    end
                end
                Hilbert = zeros(size(LFP_certain));
                for j = 1:no
                    Hilbert(j,:) = hilbert(LFP_certain(j,:));
                end
                for tt = start:nonoverlap:steps-2*bin
                    period = tt:(tt+bin-1);
                    ind = find(ismember(period,pt));
                    if ~isempty(ind)
                        t = period(ind(1));
                    else
                        continue
                    end
                    cs(1) = px(pt == t);
                    cs(2) = py(pt == t);
                    send = find(10*Distance_xy(LFP_centre_x,LFP_centre_y,cs(1),cs(2),2*hw+1) < 50);
                    if isempty(send)
                        continue
                    else
                        d = lattice_nD_find_dist(Lattice2,1.5,send);
                    end
                    receive = find(d>0)';
                    for k = 1:length(receive)
                        nte = zeros(1,bin);
                        ntesh = zeros(1,bin);
                        X = repmat(angle(Hilbert(send,t:(t+bin-1)))+pi,10,1);
                        Y = repmat(angle(Hilbert(receive(k),t:(t+bin-1)))+pi,10,1);
                        x = X(1,:);
                        y = Y(1,:);
                        Xsh = repmat(x(randperm(bin)),1,1);
                        Ysh = repmat(y(randperm(bin)),1,1);
                        for j = -(1:bin)
                            Opt.tauy = j;
                            try
                                [NTE] = transferentropy(X',Y',Opt,'NTE');
                                NTE = NTE(NTE ~= Inf);
                                NTE = NTE(NTE ~= -Inf);
                                [NTEsh] = transferentropy(Xsh',Ysh',Opt,'NTE');
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
                METHOD = [METHOD method];
                DC = [DC dc];
        end
    end
    save(['NTE5_100-',sprintf('%04g',loop_num),'.mat'],'METHOD','DC')
end
toc;
end