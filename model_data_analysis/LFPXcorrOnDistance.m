function [D,NXC,err2] = LFPXcorrOnDistance(R)
% xcorrelation of LFP on distance(normalized to power)
% tic;
% dir_strut = dir('*RYG.mat');
% num_files = length(dir_strut);
% files = cell(1,num_files);
% for id_out = 1:num_files
%     files{id_out} = dir_strut(id_out).name;
% end
[Lattice2, ~] = lattice_nD(2,1.5); % unit 150 um % 1.5
fs = 10000; % Hz
% for i = 1:num_files
%     fprintf('Loading RYG.mat file %s...\n', files{i});
%     R = load(files{i});
%     pairs = nchoosek(1:size(R.LFP.LFP_gamma,1),2);
%     for j = 1:length(pairs)
%         m = pairs(j,1);
%         n = pairs(j,2);
%         d(j) = Distance_xy(Lattice2(m,1),Lattice2(m,2),Lattice2(n,1),Lattice2(n,2),4);
%         xc = max(xcorr(R.LFP.LFP_gamma(m,:),xcorr(R.LFP.LFP_gamma(n,:))));
%         [pxxm,~] = pwelch(R.LFP.LFP_gamma(m,:),5000,3000,5000,fs); % fm 0:2:5000
%         [pxxn,~] = pwelch(R.LFP.LFP_gamma(n,:),5000,3000,5000,fs);
% %         pxxm = 10*log10(pxxm);
% %         pxxn = 10*log10(pxxn);
%         p = sum(pxxm(16:41).*pxxn(16:41)); % 30-80 Hz
%         nxc(j) = xc/p;
%     end
% %     plot(d,nxc,'o')
%
%     D = unique(d);
%     for j = 1:length(D)
%         xc = nxc(d == D(j));
%         err2(j) = std(xc);
%         NXC(j) = mean(xc);
%     end
%     errorbar(150*D,1e-8*NXC,1e-8*err2)
%     xlim([100 450])
%     xlabel('Distance(um)')
%     ylabel('Normalized X-correlation(a.u.)')
% end
dd = [];
nxcnxc = [];
% for i = 1% :2 % num_files
%     fprintf('Loading RYG.mat file %s...\n', files{i});
%     R = load(files{i});
    pairs = nchoosek(1:size(R.LFP.LFP_gamma,1),2);
    for j = 1:length(pairs)
        m = pairs(j,1);
        n = pairs(j,2);
        d(j) = Distance_xy(Lattice2(m,1),Lattice2(m,2),Lattice2(n,1),Lattice2(n,2),4); % 4 or 10
        xc = max(xcorr(R.LFP.LFP_gamma(m,:),xcorr(R.LFP.LFP_gamma(n,:))));
        [pxxm,~] = pwelch(R.LFP.LFP_gamma(m,:),5000,3000,5000,fs); % fm 0:2:5000
        [pxxn,~] = pwelch(R.LFP.LFP_gamma(n,:),5000,3000,5000,fs);
        p = sum(pxxm(16:41).*pxxn(16:41)); % 30-80 Hz
        nxc(j) = xc/p;
    end
    %     plot(d,nxc,'o')
    dd = [dd d];
    nxcnxc = [nxcnxc nxc];
% end
D = unique(dd);
for j = 1:length(D)
    xc = nxcnxc(dd == D(j));
    err2(j) = std(xc);
    NXC(j) = mean(xc);
end
% subplot(1,2,1)
errorbar(150*D,NXC,err2) % E16~150*D; E100~60*D
% xlim([100 450])
% ylim([10 14])
xlabel('Distance(um)')
ylabel('Normalized X-correlation(a.u.)')
text(-0.2,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
savefig('LFPXCorrOnDistance.fig')
% toc;
end