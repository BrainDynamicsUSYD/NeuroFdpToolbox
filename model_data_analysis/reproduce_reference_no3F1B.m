% reproduce reference NO.3 Fig1B/Fig2E
% must run reproduce_reference_no3F1A to get new R in advance
% IEI = R.LFP.Gamma_IEI;
% amp = R.LFP.Gamma_amp;

mat = zeros(42); % 46 ; 42
mat = flipud(mat);
imagesc(mat)
oldcmap = colormap(gray);
colormap( flipud(oldcmap) );
colorbar
% set(gca,'XTickLabel',[10:10:90])
% set(gca,'YTickLabel',flip([5:5:45],2))
% xlabel('LFP Amplitude(a.u.)');
% ylabel('IEI(ms)');
set(gca,'XTick',[0:6:42])
set(gca,'YTick',[0:6:42])
set(gca,'XTickLabel',[0:0.018:0.126])
set(gca,'YTickLabel',flip([0:0.162:1.134],2))
xlabel('g_E(us)');
ylabel('g_I(us)');
% [r,p] = corrcoef(amp',IEI');
[r,p] = corrcoef(G_E',G_I');
if p(2) > 10^(-4)
    pvalue = [', p = ',num2str(p(2),'%.2f')];
else
    pvalue = [', p < 10^{-4}'];
end
title(['r = ',num2str(r(2),'%.2f'),pvalue])