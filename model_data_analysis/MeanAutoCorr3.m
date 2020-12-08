function MeanAutoCorr3(varargin)
% long range correlation running in cluster, adjust from MeanAutoCorr2.m
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

% Loop number for PBS array job
loop_num = 0;

% no = size(R.LFP.LFP_gamma,1);
scales = R.LFP.wavelet.scales;
% nF = length(scales);
for window = 1e3*[150:50:1e3] % ms [5e3:5e3:1e6]
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
        % decrease the sampling rate
        LFP = mean(vec2mat(LFP,10),2)';
        for k = 10 %:nF
            coeffs_tmp = abs(cwt(LFP,scales(k),'cmor1.5-1'));
            toc;
            [env,~] = envelope(coeffs_tmp,12,'peak');            
            F = DFA(abs(env),window);
%             A = polyfit(log(dt*window),log(F),1);
%             fitline = polyval(A,log(dt*window));
%             Alpha1 = A(1);
            save(sprintf('Elec1Scales10Window%d',1e-3*window),'F','window')
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
