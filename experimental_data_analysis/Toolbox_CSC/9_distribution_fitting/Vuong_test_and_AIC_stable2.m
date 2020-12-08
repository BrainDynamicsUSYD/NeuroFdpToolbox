function [R,reject_decision,p,aic,aic_power_law] = Vuong_test_and_AIC_stable2(x)
% Vuong test of ratio of log-likelihood
% and Akaike information criterion (AIC)
% Guozhang Chen, 7/26/2018
% this function needs the toolbox of Clauset 2009
 
% fit power law
pd = fitdist(x','stable');
L_log = -sqrt(negloglik(pd));
% visualize cdf and fitting
% h = plplot(x, xmin, alpha);
% calculate power law's log-likelihoods
% L_log_power_law = -alpha*log(x(x>=xmin)) - log(zvec(find(vec<=alpha,1,'last')) - 1/n*sum((1:xmin-1).^-alpha));
% L_log_power_law = -alpha*log(x(x>=xmin));
aic_power_law = aicbic(L_log,2);
% fit other distribution and get ratio of log-likelihoods and p-value from
% Vuong test
distname = {'Normal'};
 
for ii = 1:length(distname)
    [R(ii),reject_decision(ii),p(ii),aic(ii)] = get_log_likelihoods_ratio_and_AIC(x,L_log,char(distname{ii}));
end
end
 
function [R,reject_decision,p,aic] = get_log_likelihoods_ratio_and_AIC(x,L_log_power_law,distname)
if size(x,2) > size(x,1)
    x = x';%X must be a numeric column vector.
end
pd = fitdist(x,distname);
% y = pdf(pd,x);
% figure
% plot(sort(x),pdf(pd,sort(x)))
% hold on
% histogram(x,'normalization','probability')
% set(gca,'xscale','log')
% title(distname)
% xlabel('x')
% ylabel('Probability')
L_log = -negloglik(pd);
% log likehood ratio
R = L_log_power_law - L_log;
% % Vuong test
% sigma_squre = mean((L_log_power_law - L_log - (mean(L_log_power_law) - mean(L_log))).^2);
% p = erfc(abs(R)/sqrt(2*length(x)*sigma_squre));

switch distname
    case 'Lognormal'
        numParam = 2;        
    case 'Normal'
        numParam = 2;
    case 'Poisson'
        numParam = 1;
    case 'Exponential'
        numParam = 1;
    case 'Gamma'
        numParam = 2;
end
% built in method for
[reject_decision,p] = lratiotest(L_log_power_law,L_log,numParam);
% AIC
aic = aicbic(L_log,numParam);
end