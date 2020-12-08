function PowerStateAssociation(R,P)
% load RYG.mat and pattern.mat in advance.

GammaLFP = R.LFP.LFP_gamma;
[no,tot] = size(GammaLFP);
Pt = ([1:tot-1] + [2:tot])/2;
for i = 1:no
    phase = angle(hilbert(GammaLFP(i,:)));
    PInd = find(diff(phase) < 0);
end
end