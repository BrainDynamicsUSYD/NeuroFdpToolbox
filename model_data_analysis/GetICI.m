function R = GetICI(R)
% inter-cycle-interval from LFP oscillation peak to peak,
% peak with consecutive next ICI

[no, steps] = size(R.LFP.LFP_broad);

fprintf('\t Getting IEI distribution...\n');
GammaPeak = [];
ICI = [];
flag = 0;

try
    LFP_gamma = R.LFP.LFP_gamma;
catch
    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
    % Butterworth filter
    order = 4; % 4th order
    lowFreq_br = 30;
    hiFreq_br = 80;
    Wn = [lowFreq_br hiFreq_br]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    LFP_gamma = zeros(no, steps);
    flag = 1;
end

for i = 1:no
    if flag == 1
        LFP_gamma(i,:) = filter(b,a,R.LFP.LFP{1}(i,:)); % R.LFP.LFP{1}
    end
    
    [p1,l1] = findpeaks(LFP_gamma(i,:));
    GammaPeak = [GammaPeak p1(1:end-1)]; % unit:a.u.
    ICI = [ICI 0.1*(l1(2:end) - l1(1:end-1))]; % ms    
end
% LFP_Gamma = nanmean(LFP_gamma,1);

% record results
if flag == 1
    R.LFP.LFP_gamma = LFP_gamma;
end

R.LFP.GammaICI = ICI;
R.LFP.GammaLFPPeak = GammaPeak;
end