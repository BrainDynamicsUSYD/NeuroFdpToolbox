% This function calculates the correlation between channels
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
sigGammaTemp = reshape(sigGamma,100,[]) ;
    no = size(sigGammaTemp,1);
    mat = zeros(sqrt(no));
    ref = randperm(no,1);
    timeRange = fix(1*fsTemporal)+1 : fix(3*fsTemporal) ;
    S1 = sigGammaTemp(ref,timeRange);
    for j = 1:no
        S2 = sigGammaTemp(j,timeRange);
        Corr = corrcoef(S1,S2);
        mat(j) = Corr(2);
    end
    [row,col] = ind2sub(size(mat),ref);
    imagesc(mat)
    colorbar
    title('Gamma LFP Correlation Map')
    hold on;
    plot(col,row,'k*', 'MarkerSize', 10)
