% test the influence of firing rate to calculation of transfer entropy
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 4;
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.tauy = -50;
fr = fr(1:18);% [0.1:0.1:1 1.2:0.2:3 3.5:0.5:7 8:35];
ntesh = zeros(size(fr));
for j = 1:length(fr)
    tic;
    x = rand(1,1e5) < fr(j)/1e4;
    y = rand(1,1e5) < fr(j)/1e4;
    for repeat = 1:60
        for t = 5e3:2e2:1e5
            try
                [NTEsh] = transferentropy(repmat(x((t - 499):t)',1,10),repmat(y((t - 499):t)',1,10),Opt,'NTEsh');
                NTEsh = NTEsh(NTEsh ~= Inf);
                NTEsh = NTEsh(NTEsh ~= -Inf);
                if ~isnan(nanmean(NTEsh))
                    ntesh(j) = ntesh(j) + nanmean(NTEsh);
                end
            catch
            end
        end
    end
    toc;
end