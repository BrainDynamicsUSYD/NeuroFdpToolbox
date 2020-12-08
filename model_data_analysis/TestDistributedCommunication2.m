function TestDistributedCommunication2(varargin)
tic;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

dir_strut2 = dir('pattern*.mat');
files2 = cell(1,num_files);
for id_out = 1:num_files
    files2{id_out} = dir_strut2(id_out).name;
end

%%% basic setting %%%
hw = 31;
spots = 7;
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 4;
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 10; % number of calculating trails
Opt.testMode = 1;
nX1 = [];
nX2 = [];
nY1 = [];
nY2 = [];
dmin = 0;
[Lattice1, ~] = lattice_nD(2, hw);

%%% Loop number for PBS array job %%%
loop_num = 0;

%%%%%%%%%%%%%%% variables %%%%%%%%%%%%%%%%%%
bin = 300; % 0.1ms
start = 10006; % R.grid.t_mid = 26:10:99966

% mode = 1; % 1:calculate spiking trains; 2:LFP phase
% method = 1; % 1:random select sites; 2:pattern site as sender

%%% mode 1:spikes %%%

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
    [~, steps] = size(R.LFP.LFP_broad);
    pt = 10*P.ts;
    px = P.x;
    py = P.y;
    M = ones(1,1)*[11 12];
    for method = (M(:))' % [11 12]
        
        % set seed
        date_now = datestr(now,'yyyymmddHHMM-SSFFF');
        scan_temp = textscan(date_now,'%s','Delimiter','-');
        rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
        alg='twister';
        rng(rand_seed,alg); %  Note that this effect is global!
        
        switch method
            case 11
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
                for tt = start:bin:steps-bin
                    period = tt:(tt+bin-1);
                    ind = find(ismember(period,pt));
                    if ~isempty(ind)
                        t = period(ind(1));
                    else
                        continue
                    end
                    X = R.spike_hist{1}(FindNeurons(Lattice1,hw,10,send)',(t-bin+1):t);
                    for r = 1:spots-1
                        Y = R.spike_hist{1}(FindNeurons(Lattice1,hw,10,receive(r))',(t-bin+1):t);
                        for j = -(1:bin)
                            Opt.tauy = j;
                            try
                                [NTE] = transferentropy(X',Y',Opt,'NTE');
                                NTE = NTE(NTE ~= Inf);
                                NTE = NTE(NTE ~= -Inf);
                                if nanmean(NTE) > 0
                                    nX1 = [nX1 sum(X(:))];
                                    nY1 = [nY1 sum(Y(:))];
                                    break
                                end
                            catch
                            end
                        end
                    end
                end
            case 12
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
                for tt = start:bin:steps-bin
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
                    X = R.spike_hist{1}(FindNeurons(Lattice1,hw,10,send)',(t-bin+1):t);
                    for r = 1:spots-1
                        Y = R.spike_hist{1}(FindNeurons(Lattice1,hw,10,receive(r))',(t-bin+1):t);
                        for j = -(1:bin)
                            Opt.tauy = j;
                            try
                                [NTE] = transferentropy(X',Y',Opt,'NTE');
                                NTE = NTE(NTE ~= Inf);
                                NTE = NTE(NTE ~= -Inf);
                                if nanmean(NTE) > 0
                                    nX2 = [nX2 sum(X(:))];
                                    nY2 = [nY2 sum(Y(:))];
                                    break
                                end
                            catch
                            end
                        end
                    end
                end
        end
    end
    save(['teen',sprintf('%04g',loop_num),'.mat'],'nX1','nY1','nX2','nY2')
end
toc;
end