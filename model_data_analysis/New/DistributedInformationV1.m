function DistributedInformationV1(varargin)
% based on signal collect from LFPAmpPattern.m function
% adapt from DistributedCommincation2.m function
% select the most 6 frequent electrodes in pattern as the parallel signal
% communication between purely patterns
% DC Index in certain interval
tic;
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
dir_strut2 = dir('3DBurst0*minTime15SR1000P95.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end

%%% basic setting %%%
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 5; % number of trials - 1
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 6; % number of calculating trails
Opt.testMode = 1;

%%% Loop number for PBS array job %%%
loop_num = 0;

bin = 10;
dt = 0.1;
for i = 1:num_files
    for interval = [2.1e4:1e3:3.4e4]/bin
        for gap = [1e4]/bin
            
            % For PBS array job
            loop_num = loop_num + 1;
            if nargin ~= 0
                PBS_ARRAYID = varargin{1};
                if loop_num ~=  PBS_ARRAYID
                    continue;
                end
            end
            
            fprintf('Loading RYG.mat file %s...\n', files{i});
            R = load(files{i});
            fprintf('Loading RYG.mat file %s...\n', files2{i-9});
            P = load(files2{i-9});
            [no,steps] = size(R.LFP.LFP{1});
            LFP_certain = zeros(no,steps/bin);
            LFP_suffle = zeros(no,steps/bin);
            DC = cell(1);
            
            % set seed
            date_now = datestr(now,'yyyymmddHHMM-SSFFF');
            scan_temp = textscan(date_now,'%s','Delimiter','-');
            rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
            alg='twister';
            rng(rand_seed,alg); %  Note that this effect is global!
            
            fs = 1/(dt*bin*1e-3); % sampling frequency (Hz)
            
            % Butterworth filter
            order = 4; % 4th order
            lowFreq_br = 5;
            hiFreq_br = 100;
            Wn = [lowFreq_br hiFreq_br]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
            for n = 1:no
                LFP_certain(n,:) = filter(b,a,R.LFP.LFP{1}(n,1:bin:end)); % R.LFP.LFP{1}
                LFP_suffle(n,:) = ifft(fft(LFP_certain(n,:)).*exp(1j*pi*rand(1,steps/bin)),'symmetric');
            end
            Hilbert = zeros(size(LFP_certain));
            HilbertS = zeros(size(LFP_certain));
            for n = 1:no
                Hilbert(n,:) = hilbert(LFP_certain(n,:));
                HilbertS(n,:) = hilbert(LFP_suffle(n,:));
            end
            np = length(P.ts);
            for k = 1:np
                send = sort(P.Send{k},'descend');
                P.Send{k} = send(1:6);
            end
            p = 1:gap:(steps/bin-interval+1);
            for per = 1:length(p)
                cand = find(P.ts>=p(per) & P.ts+P.Duration<=p(per)+interval-1);
                dc = [];
                for k = 1:length(cand)-1
                    l = P.Duration(cand(k));
                    nte = zeros(1,l);
                    ntesh = zeros(1,l);
                    % LFP phase as candidate signal
                    %             X = repmat(angle(Hilbert(send,t:(t+bin-1)))+pi,4,1);
                    %             Y = repmat(angle(Hilbert(receive(k),t:(t+bin-1)))+pi,4,1);
                    %             x = X(1,:);
                    %             y = Y(1,:);
                    %             Xsh = repmat(x(randperm(bin)),4,1);
                    %             Ysh = repmat(y(randperm(bin)),4,1);
                    % LFP amplitude as candidate signal
                    X = abs(Hilbert(P.Send{cand(k)},P.ts(cand(k)):(P.ts(cand(k))+l-1)));
                    Y = abs(Hilbert(P.Send{cand(k+1)},P.ts(cand(k)):(P.ts(cand(k))+l-1)));
                    Xsh = abs(HilbertS(P.Send{cand(k)},P.ts(cand(k)):(P.ts(cand(k))+l-1)));
                    Ysh = abs(HilbertS(P.Send{cand(k+1)},P.ts(cand(k)):(P.ts(cand(k))+l-1)));
                    for tau = -(1:l)
                        Opt.tauy = tau;
                        try
                            [NTE] = transferentropy(X',Y',Opt,'NTE');
                            NTE = NTE(NTE ~= Inf);
                            NTE = NTE(NTE ~= -Inf);
                            [NTEsh] = transferentropy(Xsh',Ysh',Opt,'TE','NTE');
                            NTEsh = NTEsh(NTEsh ~= Inf);
                            NTEsh = NTEsh(NTEsh ~= -Inf);
                            ntesh(abs(tau)) = nanmean(NTEsh);
                            if nanmean(NTE) > 0
                                nte(abs(tau)) = nanmean(NTE);
                            end
                        catch
                        end
                    end
                    if max(nte) > max([max(ntesh),0])
                        dc = [dc max(nte)];
                    end
                end
                DC{per} = dc;
            end
            save(['DCAmplitudeSAmplitudePmin15msInt',num2str(interval),...
                'msGap',num2str(gap),'ms-',sprintf('%04g',i),'.mat'],'DC')
        end
        toc;
    end
end
end