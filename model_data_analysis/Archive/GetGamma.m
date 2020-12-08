function [R] = GetGamma(R)
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.
% adjust from get_GammaV1 function
fprintf('\t Getting Gamma band filtered signal...\n');
dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)

LFP = R.LFP.LFP{1};
[no,~] = size(LFP);

% Butterworth filter
order = 4; % 4th order
lowFreq_br = 1; % broad band (1-1000 Hz)
hiFreq_br = 1000;
Wn = [lowFreq_br hiFreq_br]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
LFP_broad = zeros(size(LFP));
for i = 1:no
    LFP_broad(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
end

% Butterworth filter
order = 4; % 4th order
R.LFP.lowFreq = 30; % gamma band
R.LFP.hiFreq = 80;
Wn = [R.LFP.lowFreq R.LFP.hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
LFP_gamma_hilbert_abs = zeros(size(LFP));
gaus_width = 12.5; % ms
[ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );
LFP_gamma = zeros(size(LFP)); 
for i = 1:no
    LFP_gamma(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
%     LFP_gamma_hilbert_abs(i,:) = abs(hilbert(LFP_gamma(i,:)));
%     % hilbert transformation & gaussian smoothing
    LFP_gamma_hilbert_abs(i,:) = conv(abs(hilbert(LFP_gamma(i,:))), Kernel,'same'); %#ok<AGROW>
end

% Output results
R.LFP.LFP_broad = LFP_broad;
R.LFP.LFP_gamma = LFP_gamma; % LFP_gamma
R.LFP.LFP_gamma_hilbert_abs = LFP_gamma_hilbert_abs;
R.LFP.gauss_width = gaus_width;
clear LFP_broad LFP_gamma LFP_gamma_hilbert;

freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
Fs = 1000/dt;
fc = centfrq('cmor1.5-1');
scalerange = fc./(freqrange*(1/Fs));
scales = scalerange(end):0.5:scalerange(1);
pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies

R.LFP.wavelet.pseudoFreq = pseudoFreq;
R.LFP.wavelet.scales = scales;
fprintf('\t Finishing get_Gamma...\n');
end