% long range correlation
% average autocorrelation
% DFA
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
dt = 0.1; % ms
max_lag = 1e6; % ms
AC = [];
window = [5e4:5e4:3e5]; % 0.1 ms
for i = 3 % 1%:num_files % [1 12 17 18]
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i},'LFP');
    disp('done.\n');
    no = size(R.LFP.LFP_gamma,1);
    scales = R.LFP.wavelet.scales;
    nF = length(scales);
    for j = 1%:no
        LFP = R.LFP.LFP_gamma(j,:);
        for k = 1 % 50:50:nF %:nF
            coeffs_tmp = abs(cwt(LFP,scales(k),'cmor1.5-1'));
            [env,~] = envelope(coeffs_tmp,30,'peak');
            
            F = zeros(size(window));
            for tau = 1:length(window)
                F(tau) = DFA(abs(env),window(tau));
            end
            figure
            loglog(dt*window,F,'o');
%             plot(log(dt*window),log(F));
            xlabel('Window Size(ms)')
            ylabel('F(tau)')
            A = polyfit(log(dt*window),log(F),1);
            hold on
            fitline = polyval(A,log(dt*window));
            loglog(dt*window,exp(fitline))
            Alpha1 = A(1);
            disp(Alpha1)
%                     plot(1:1e5,LFP,1:1e5,LFPenvelope)
%                         [ac,lags] = autocorr(env, round(max_lag/dt) );
%                         AC = [AC;ac];
        end
    end
end
% Autoc = mean(AC,1);
% % Autoc = tsmovavg(Autoc,'s',1e4,2);
% plot(lags*dt*1e-3,Autoc,'k')
% xlabel('Time(s)')
% ylabel('Correlation')

plot(log(dt*1e-3*window),log(F));
xlabel('Window Size(s)')
ylabel('F(tau)')
A = polyfit(log(window),log(F),1);
Alpha1 = A(1);