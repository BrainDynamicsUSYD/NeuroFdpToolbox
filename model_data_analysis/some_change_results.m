function some_change_results(R)
figure(1)
subplot(8,1,[1:7])
prn = cell2mat( R.LFP.wavelet.peak.prn(~cellfun(@isempty,R.LFP.wavelet.peak.prn)) );
rp_raw_amp = cell2mat( R.LFP.wavelet.peak.rp_raw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.rp_raw_amp)) );
rp_freq = cell2mat( R.LFP.wavelet.peak.rp_freq(~cellfun(@isempty,R.LFP.wavelet.peak.rp_freq)) );
X = [rp_raw_amp(~isnan(prn))',prn(~isnan(prn))'];
[N, C] = hist3(X,[20 20]);
imagesc(C{1}, C{2},N);
set(gca, 'YDir','normal')
colorbar;
xlabel('Ripple magnitude')
ylabel('Post-ripple negativtiy (mininum LFP)')
Post_ripple_negativity = corrcoef(rp_raw_amp(~isnan(prn)),prn(~isnan(prn)));
Post_ripple_negativity = Post_ripple_negativity(1,2);
title(['r = ' num2str(Post_ripple_negativity)])
com = subplot(8,1,8);
set(com,'Visible','off'); 
comments = sprintf(['ripple frequency:%0.5g ± %0.5g Hz\nSWRs frequency:'...
    '%0.5g ± %0.5g Hz'],mean(rp_freq),std(rp_freq),mean(R.LFP.ripple_event.Hz),std(R.LFP.ripple_event.Hz));
text(0,-0.5,comments);
end