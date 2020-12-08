function MaximalFrequency
% find the instantenous maximal frequency in wavelet time series, i.e. the
% frequency which has maximum CData in wavelet analysis
tic
dir_strut = dir('*_out_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
for n = 1 % :num_files % 1:num_files
    fprintf('Loading out_RYG.mat file %s...\n', files{n});
    R = load(files{n});
    dt = R.dt;
    step_tot = R.step_tot;
    seg_size = 2e4; % 2s
    seg_num = ceil(step_tot/seg_size);
    [no, ~] = size(R.LFP.LFP_gamma);
    scales = R.LFP.wavelet.scales;
    pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies
    maxFre = [];
    for i = 1:no
        maxF = [];
        for seg = 1:seg_num
            seg_ind = get_seg(step_tot, seg_size, seg);
%             t = seg_ind*dt*1e-3; % second           
            x_tmp = R.LFP.LFP_gamma(i,seg_ind);
            coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
            CData = transpose(coeffs_tmp); % coeffs_tmp'
            [~,I] = max(CData);
            maxF = [maxF pseudoFreq(I)];
        end
        maxFre = [maxFre;maxF];
    end
    mf = maxFre(logical(R.LFP.GammaBurstEvent.is_burst));
end
histogram(mf)
toc
end
