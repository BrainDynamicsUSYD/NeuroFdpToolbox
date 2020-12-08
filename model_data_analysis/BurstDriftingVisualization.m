function BurstDriftingVisulization(R,no_given,seg_given)

dt = R.dt;
step_tot = R.step_tot;
seg_size = 2e4; % 2s
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
for i = 1 % nos
    if nargin == 3
        seg_num = 1;
    end
    for seg = 4 % 1:seg_num
        seg_ind = get_seg(step_tot, seg_size, seg);
        if nargin == 3
            seg_ind = seg_given;
        end
        t = (seg_ind-1)*dt*1e-3; % second
        fsTemporal = 1/(dt*1e-3); % sampling frequency (Hz)
        
        x_tmp = R.LFP.LFP_gamma(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
        %         fc = centfrq('cmor1-5');
        %         scalerange = fc./([30 80]*(1/fs));
        %         scales = scalerange(end):0.5:scalerange(1);
        %         pseudoFreq = scal2frq(scales,'cmor1-5',1/fs);
        %         coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1-5'))'; % this had narrow width of scales
        CData = transpose(coeffs_tmp); % coeffs_tmp'
%         zscoreCoef = zscore(CData,1);        
%         subplot(2,1,1)
%         uimagesc(t,pseudoFreq(end:-1:1),zscoreCoef(end:-1:1,:));
% %         uimagesc(t,pseudoFreq(end:-1:1),CData(end:-1:1,:));
%         colorbar('off');
%         xlim([2*(seg-1) 2*seg]);
%         ylim(frequency);
%         xlabel('Time(s)');
%         ylabel('Frequency(Hz)');
%         set(gca,'YDir','normal')
        
%         subplot(2,1,2)
        zscoreCoef = CData;
        Y = prctile(zscoreCoef(:),95) ;
        binaryImage = zeros(size(zscoreCoef)) ;
        binaryImage(zscoreCoef>=Y) = 1 ;
        CC = bwconncomp(binaryImage(:,1:end)) ;
        validRegionIdx = 0;
        count = 1;
        for iRegion = 1:size(CC.PixelIdxList,2)
            if(size(CC.PixelIdxList{iRegion},1)>200)
                validRegionIdx(count) = iRegion ;
                count = count + 1 ;
            end
        end
        temp = zeros(size(binaryImage(:,1:2*fsTemporal)));
        for iRegion = 1:length(validRegionIdx)
            temp(CC.PixelIdxList{validRegionIdx(iRegion)}) = 1 ;
        end
        uimagesc(t,pseudoFreq(end:-1:1),temp(end:-1:1,:))
        % 2000: 1Hz 400: 5Hz 60:30Hz 20:80Hz
        set(gca,'YDir','normal')
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        text(-0.1,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
        
%         next = input('\t Next figure?');
%         delete(gcf);
    end
    toc;
end
end