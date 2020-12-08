function MeanAutoCorr2(varargin)
% long range correlation running in cluster, adjust from MeanAutoCorr.m
% average autocorrelation
% DFA
R = load('0003-201804120739-17816_in_1523482960796_out_RYG.mat','LFP');
% dir_strut = dir('*RYG.mat');
% num_files = length(dir_strut);
% files = cell(1,num_files);
% for i = 1:num_files
%     files{i} = dir_strut(i).name;
% end
% R = load(files{3},'LFP');
disp('Loading done.\n');
dt = 0.1; % ms

% Loop number for PBS array job
loop_num = 0;

% no = size(R.LFP.LFP_gamma,1);
scales = R.LFP.wavelet.scales;
% nF = length(scales);
for window = [5e4:5e4:1e6] % 0.1 ms
    loop_num = loop_num + 1;    
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    tic;
    for j = 1 % :no
        LFP = R.LFP.LFP_gamma(j,:);
        for k = 10 %:nF
            coeffs_tmp = abs(cwt(LFP,scales(k),'cmor1.5-1'));
            toc;
            [env,~] = envelope(coeffs_tmp,30,'peak');            
            F = DFA(abs(env),window);
%             A = polyfit(log(dt*window),log(F),1);
%             fitline = polyval(A,log(dt*window));
%             Alpha1 = A(1);
            save(['Elec1Scales10Window%d',1e-4*window],'F')
            toc;
            %             loglog(dt*window,F,'o');
            %             xlabel('Window Size(ms)')
            %             ylabel('F(tau)')
            %             hold on
            %             loglog(dt*window,exp(fitline))
            %             disp(Alpha1)
        end
    end
end
end
