function WaveletImagescOnPhase
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
phase = -pi:2*pi/50:pi;
for i = 1 % :num_files
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    % Butterworth filter
    order = 4;
    lowFreq = 4;
    hiFreq = 10;
    fs = 1e4;
    Wn = [lowFreq hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
    seg_size = 1e4; % 1s
    seg_num = ceil(R.step_tot/seg_size);
    for j = 1:size(R.LFP.LFP{1},1)
        mat = zeros(length(R.LFP.wavelet.scales),50,seg_num-1);
        for seg = 2:seg_num
            seg_ind = get_seg(R.step_tot, seg_size, seg);
            LFP_theta = angle(hilbert(filter(b,a,R.LFP.LFP{1}(j,seg_ind))));
            x_tmp = R.LFP.LFP_gamma(j,seg_ind);
            CData = abs(cwt(x_tmp,R.LFP.wavelet.scales,'cmor1.5-1'));                  
            for k = 1:50
                ind = find(LFP_theta>=phase(k) & LFP_theta<=phase(k+1));
                mat(:,k,seg-1) = mean(CData(:,ind),2);
            end
        end
        mat = mean(mat,3);
        pv = (phase(1:end-1) + phase(2:end))/2;
        subplot(4,4,j)
        uimagesc(pv,R.LFP.wavelet.pseudoFreq(end:-1:1),mat(end:-1:1,:));
        set(gca,'YDir','normal')
        xlabel('Theta Phase(rad)')
        ylabel('Frequency(Hz)')
%         next = input('\t Next figure?');
%         delete(gcf);
    end
end
end