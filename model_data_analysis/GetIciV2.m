function [R] = GetIciV2(R)
% inter-cycle-interval from LFP oscillation peak to peak

dt = R.dt;
% [no, steps] = size(R.LFP.LFP_gamma); % R.LFP.LFP_broad
[no, steps] = size(R.LFP.LFP_broad);

fprintf('\t Getting IEI distribution...\n');
IEI_lumped = [];
Ind = [];
IEI_lumped2 = [];
Ind2 = [];
LFP_peak = [];
LFP_valley = [];
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
    
    LFP_peak_t = []; % time steps
    LFP_valley_t = [];
    for j = 2:(steps - 1)
        if LFP_gamma(i,j) > LFP_gamma(i,j - 1) && LFP_gamma(i,j) > LFP_gamma(i,j + 1)
            peak = j;
            LFP_peak_t = [LFP_peak_t peak];
        end
        if LFP_gamma(i,j) < LFP_gamma(i,j - 1) && LFP_gamma(i,j) < LFP_gamma(i,j + 1)
            valley = j;
            LFP_valley_t = [LFP_valley_t valley];
        end
    end
    LFP_peak = [LFP_peak LFP_peak_t];
    LFP_valley = [LFP_valley LFP_valley_t];
    % get IEI
    if length(LFP_peak_t) >= 2
        IEI_temp = (LFP_peak_t(2:end) - LFP_peak_t(1:end-1))*dt; % in ms
        ind = length(IEI_temp);
        IEI_lumped = [IEI_lumped IEI_temp];
        Ind = [Ind ind];
    end
    if length(LFP_valley_t) >= 2
        temp = (LFP_valley_t(2:end) - LFP_valley_t(1:end-1))*dt; % in ms
        ind2 = length(temp);
        IEI_lumped2 = [IEI_lumped2 temp];
        Ind2 = [Ind2 ind2];
    end
end
% LFP_Gamma = nanmean(LFP_gamma,1);

% record results
if flag == 1
    R.LFP.LFP_gamma = LFP_gamma;
end

R.LFP.GammaPeakIei = IEI_lumped;
R.LFP.GammaPeakIeiInd = Ind;
R.LFP.GammaValleyIei = IEI_lumped2;
R.LFP.GammaValleyIeiInd = Ind2;
R.LFP.GammaLfpPeak = LFP_peak;
R.LFP.GammaLfpValley = LFP_valley;
end