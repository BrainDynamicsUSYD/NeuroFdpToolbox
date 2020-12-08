function FrequencyDurationDistribution3(varargin)
% Frequency-Duration Distribution
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
% REFERENCE:Local generation and propagation of ripples along the septotemporal axis of the hippocampus
tic
R = load('0001-201702091535-53644_in_1486619228380_out_RYG.mat');

dt = R.dt;
step_tot = R.step_tot;
seg_size = 1e4; % 1s
overlap = 9e3;
seg_num = ceil((step_tot - seg_size)/(seg_size - overlap) + 1);
[R] = get_IEI(R);
[no, ~] = size(R.LFP.LFP_gamma);

freqrange = [5 100];
Fs = 1000/dt;
fc = centfrq('cmor1.5-1');
scalerange = fc./(freqrange*(1/Fs));
scales = scalerange(end):0.5:scalerange(1);
pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies from 100 to 5
fre_dur = cell(1,length(pseudoFreq));
% loop_num = 0;
toc
for i = 1:no
    for seg = 6:seg_num
        %         loop_num = loop_num + 1;
        %         % For PBS array job
        %     if nargin ~= 0
        %         PBS_ARRAYID = varargin{1};
        %         if loop_num ~=  PBS_ARRAYID
        %             continue;
        %         end
        %     end
        seg_ind = ((seg-1)*(seg_size - overlap)+1):((seg-1)*(seg_size - overlap)+seg_size);
        t = (seg_ind-1)*dt*1e-3; % second
        x_tmp = R.LFP.LFP_gamma(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        min_C = min(transpose(coeffs_tmp));
        min_C = min(min_C(:));
        max_C = max(transpose(coeffs_tmp));
        max_C = max(max_C(:));
        threshold = 200;
        LFPspec = transpose(coeffs_tmp);
        LFPspec(LFPspec < threshold) = 0;
        LFPspec(LFPspec >= threshold) = 1;
%         toc
        for j = 1:length(pseudoFreq)
            if LFPspec(j,end) ~= 0
                adjust = max(find(LFPspec(j,:) == 0));
                LFPspec(j,adjust:end) = 0;
            end
            if LFPspec(j,1) ~= 0
                adjust = min(find(LFPspec(j,:) == 0));
                LFPspec(j,1:adjust) = 0;
            end
            a = LFPspec(j,:);
            b = diff(a);
            c = find(b);
            d = diff([0,c]);
            dur = dt*d(b(c) == -1 & a(c) == 1);
            fre_dur{j} = [fre_dur{j} dur];
        end
    end
    toc
end
save('FrequencyDuration1.mat','fre_dur');
toc
end