function TimeFrequencyCImage(R,no_given,seg_given)
% frequency-time color image: Morlet wavelet spectrogram of signal(LFP)
% -1 for no figure, 0 for displaying figure, 1 for saving figure
% REFERENCE:Local generation and propagation of ripples along the septotemporal axis of the hippocampus
% adjust from function plot_SWR.m

dt = R.dt;
step_tot = R.step_tot;
seg_size = 3e4; % 1s
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_gamma);
frequency = [R.LFP.lowFreq R.LFP.hiFreq];
scales = R.LFP.wavelet.scales;
pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies
nos = 1:no;
if nargin >= 2
    nos = no_given;
end

tic;
for i = nos
    if nargin == 3
        seg_num = 1;
    end
    for seg = 1:seg_num
        seg_ind = get_seg(step_tot, seg_size, seg);
        if nargin == 3
            seg_ind = seg_given;
        end
        t = (seg_ind-1)*dt*1e-3; % second
        x_tmp = R.LFP.LFP_gamma(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        CData = transpose(coeffs_tmp);
        imagesc('XData',t,'YData',pseudoFreq,'CData',CData);
        colorbar;
        xlim([3*(seg-1) 3*seg]);
        ylim(frequency);
        xlabel('Time(s)');
        ylabel('Frequency(Hz)');
        %         name = sprintf(['%04d_Gamma_%d_%d.pdf'],R.ExplVar.loop_num,i,seg);        
        next = input('\t Next figure?');
        delete(gcf);
    end
    toc;
end
end
