function DistributedInformationTransfer(varargin)


close all;
clc;
tic;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
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
c = 1;
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
v = 31.7; % um/ms = mm/s
% speed = 16.7; % v-30:0.5:v-15;

bin = 2e3; % 0.1ms
nonoverlap = 1e3; % 0.1ms

% mode = 1; % 1:calculate spiking trains; 2:LFP phase
% method = 1; % 1:random select sites; 2:pattern site as sender


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
    R = get_grid_firing_centre(R);
    [no, steps] = size(R.LFP.LFP_broad);
    start = R.grid.t_mid(1e3); % R.grid.t_mid = 26:10:99966
    
    for mode = [1 2]
        for method = [1 2]
            ntesh1 = 0; % zeros(1,180); % num_files);
            ntesh2 = 0;
            num1 = 0; % zeros(1,180); % num_files);
            num2 = 0;
            
            
            % set seed
            date_now = datestr(now,'yyyymmddHHMM-SSFFF');
            scan_temp = textscan(date_now,'%s','Delimiter','-');
            rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
            alg='twister';
            rng(rand_seed,alg); %  Note that this effect is global!
            
%             cs = round(R.grid.quick.centre);
%             I = cs(1,:);
%             J = cs(2,:);
%             csi = ~isnan(I).*(ismember(I,-31:31)).*(ismember(J,-31:31));
%             csind = [1:length(csi)].*csi;
%             csind = csind(csind > 0);
%             tc = R.grid.t_mid(csind);
%             I = I(csind);
%             J = J(csind);
%             sendc = IndexLattice(Lattice1,I,J);
%             [a,~] = hist(sendc,1:3969);
%             send = find(a == max(a));
%             send = send(1);
%             indt = find(sendc == send);
            
            
            if mode == 1
                b = randperm(length(Lattice1),spots);
                receive = b(2:spots);
%                 for t = tc(indt)
                for t = start:nonoverlap:steps
                    if method == 1
                        send = b(1);
                    else
                        cs = round(R.grid.quick.centre(:,R.grid.t_mid == t));
                        if sum(isnan(cs)) > 0 || max(abs(cs)) > hw
                            continue
                        else
                            send = IndexLattice(Lattice1,cs(1),cs(2));
                        end
                    end
                    tauy = -round(Distance(Lattice1(receive,1),Lattice1(receive,2),Lattice1(send,1),Lattice1(send,2),hw*2 + 1)*10/v*10); % 0.1ms
                    for r = 1:spots-1
                        Opt.tauy = tauy(r); % 0.1ms
                        try
                            [NTEsh] = transferentropy(R.spike_hist{1}(FindNeurons(Lattice1,hw,10,send)',(t-bin+1):t)',R.spike_hist{1}(FindNeurons(Lattice1,hw,10,receive(r))',(t-bin+1):t)',Opt,'NTEsh');
                            NTEsh = NTEsh(NTEsh ~= Inf);
                            NTEsh = NTEsh(NTEsh ~= -Inf);
                            if nanmean(NTEsh) > 0
                                ntesh1 = ntesh1 + nanmean(NTEsh);
                                num1 = num1 + 1;
                                for nt = 1:bin
                                    En(nt) = entropy(double(full(R.spike_hist{1}(FindNeurons(Lattice1,hw,10,send)',t-bin+nt)')),Opt,'HR');
                                end
                                En(En<0) = 0;
                                if nanmean(En) > 0 && nanmean(NTEsh) <= nanmean(En)
                                    ntesh2 = ntesh2 + nanmean(NTEsh)/nanmean(En);
                                    %                             nteshm = nteshm + nanmean(NTEsh)/max(En);
                                    num2 = num2 + 1;
                                end
                            end
                        catch
                        end
                    end
                end
            else % mode 2
                try
                    LFP_gamma = R.LFP.LFP_gamma;
                catch
                    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
                    % Butterworth filter
                    order = 4; % 4th order
                    lowFreq_br = 5;
                    hiFreq_br = 100;
                    Wn = [lowFreq_br hiFreq_br]/(fs/2);
                    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
                    for j = 1:no
                        LFP_gamma(j,:) = filter(b,a,R.LFP.LFP{1}(j,:)); % R.LFP.LFP{1}
                    end
                end
                for j = 1:no
                    Hilbert(j,:) = hilbert(LFP_gamma(j,:));
                end
                
%                 ds = Distance(LFP_centre_x,LFP_centre_y,Lattice1(send,1),Lattice1(send,2),2*hw+1);
%                 send = find(ds == min(ds));
%                 send = send(1);
%                 
%                 for t = tc(indt)
                for t = start:nonoverlap:steps
                    if method == 1
                        send = No;
                    else
                        cs = round(R.grid.quick.centre(:,R.grid.t_mid == t));
                        send = find(10*Distance(LFP_centre_x,LFP_centre_y,cs(1),cs(2),2*hw+1) < 50);
                    end
                    if isempty(send)
                        continue
                    else
                        d = lattice_nD_find_dist(Lattice2,1.5,send);
                    end
                    for k = find(d>0)'
                        Opt.tauy = -round(d(k)*120/v*10); % 0.1ms
                        try
                            [NTEsh] = transferentropy(repmat(angle(Hilbert(send,(t-bin+1):t))+pi,10,1)',repmat(angle(Hilbert(k,(t-bin+1):t))+pi,10,1)',Opt,'NTEsh');
                            NTEsh = NTEsh(NTEsh ~= Inf);
                            NTEsh = NTEsh(NTEsh ~= -Inf);
                            if nanmean(NTEsh) > 0
                                ntesh1 = ntesh1 + nanmean(NTEsh);
                                num1 = num1 + 1;
                            end
                        catch
                        end
                    end
                end
            end
            ntesh(c) = ntesh1;
            num(c) = num1;
            ntesh(c+1) = ntesh2;
            num(c+1) = num2;
            c = c + 2;
            toc;
        end
    end
    % (mode,method,num)
    % 1 (1,1,1); 2 (1,1,2); 3 (1,2,1); 4 (1,2,2); 5 (2,1,1); 6 (2,1,2); %% 7 (2,2,1); 8 (2,2,2)
    
    save(['teen',sprintf('%04g',loop_num),'.mat'],'ntesh','num')
    toc;
end
end


%%% good but useless commands %%%
%         [Ind,~] = find_first_spot_grid(R);
%         candidate = round(R.grid.quick.centre(:,Ind));
%         Candidate = unique(candidate','rows');
%         mark = zeros(1,length(Candidate));
%         for j = 1:length(Candidate)
%             mark(j) = sum(ismember(candidate',Candidate(j,:),'rows'));
%         end
%         [~,b] = sort(mark,'descend');
%         try
%             send = IndexLattice(Lattice,Candidate(b(1),1),Candidate(b(1),2));
%             receive = IndexLattice(Lattice,Candidate(b(2:spots),1),Candidate(b(2:spots),2));
%         catch
%             warning('No enough candidates. Using random candidates.\n');
%         end