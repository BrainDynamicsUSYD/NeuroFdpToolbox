%% Maximum likelihood estimation
%power law distribution
pf_power = @(x,alpha,a) x.^-alpha*(alpha-1)*a^(alpha-1);
[lambdaHat,lambdaCI] = mle(Duration, 'pdf',pf_power, 'start',[1.1,30], 'lowerbound',[1,10], 'upperbound', [3,40]) ;
L = sum(log(pf_power(Duration,lambdaHat(1), lambdaHat(2) ))) ;

% exponential distribution
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(Duration, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
L3 = sum(log(pf_exp(Duration,lambdaHat3))) ;

%% MLE for 1 variable (duration)

close all
figure
xmin = 30 ;
sigIn = Duration ;
sigIn = sigIn(sigIn>xmin) ;

% power distribution with a start (xmin)
pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
[lambdaHat1,lambdaCI] = mle(sigIn, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
L1 = sum(log(pf_power(sigIn,lambdaHat1 ))) ;

z = linspace(20,1000,500);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_power(x,lambdaHat1))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat1));

% MLE of power law with cut off
figure;
pf_powerC = @(x,alpha,lambda) x.^-alpha.*exp(-lambda*x)*lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
[lambdaHat2,lambdaCI] = mle(sigIn, 'pdf',pf_powerC, 'start',[1,0], 'lowerbound',[1,0], 'upperbound', [3,10]) ;
L2 = sum(log(pf_powerC(sigIn,lambdaHat2(1), lambdaHat2(2) ))) ;


histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_powerC(x,lambdaHat2(1),lambdaHat2(2)))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_powerC(x, lambdaHat2(1),lambdaHat2(2)));


% MLE of exponential
figure;
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;

histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_exp(x,lambdaHat3))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));

%%
figure
c = histcounts(Duration,linspace(20,200,20),'Normalization','pdf');
x = linspace(20,200,20) ;
x = 0.5*(x(1:end-1) + x(2:end)) ;
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat));



%% MLE of power law for scale
close all
figure
xmin = 150 ;
sigIn = patternScale ;
sigIn = sigIn(sigIn>xmin) ;

pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
[lambdaHat1,lambdaCI] = mle(sigIn, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
L1 = sum(log(pf_power(sigIn,lambdaHat1 ))) ;

z = linspace(1e1,1e5,1000);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_power(x,lambdaHat1))

figure
c = histcounts(sigIn,linspace(1e1,1e5,2000),'Normalization','pdf');
x = linspace(1e1,1e5,2000) ;
x = 0.5*(x(1:end-1) + x(2:end)) ;
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat1));

% MLE of power law with cut off
figure;
pf_powerC = @(x,alpha,lambda) x.^-alpha.*exp(-lambda*x)*lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
[lambdaHat2,lambdaCI] = mle(sigIn, 'pdf',pf_powerC, 'start',[1,0], 'lowerbound',[1,0], 'upperbound', [3,10]) ;
L2 = sum(log(pf_powerC(sigIn,lambdaHat2(1), lambdaHat2(2) ))) ;

z = linspace(1e1,1e5,1000);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_powerC(x,lambdaHat2(1),lambdaHat2(2)))

figure
c = histcounts(sigIn,linspace(1e1,1e5,1000),'Normalization','pdf');
x = linspace(1e1,1e5,1000) ;
x = 0.5*(x(1:end-1) + x(2:end)) ;
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_powerC(x, lambdaHat2(1),lambdaHat2(2)));


% MLE of power law with cut off
figure;
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;

z = linspace(1e1,1e5,1000);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_exp(x,lambdaHat3))

figure
c = histcounts(sigIn,linspace(1e1,1e5,1000),'Normalization','pdf');
x = linspace(1e1,1e5,1000) ;
x = 0.5*(x(1:end-1) + x(2:end)) ;
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));


%% MLE for inverval


close all
figure
xmin = 50 ;
sigIn = centInterval ;
sigIn = sigIn(sigIn>xmin) ;

pf_power = @(x,alpha) x.^-alpha*(alpha-1)*xmin^(alpha-1);
[lambdaHat1,lambdaCI] = mle(sigIn, 'pdf',pf_power, 'start',1, 'lowerbound',1, 'upperbound', 3) ;
L1 = sum(log(pf_power(sigIn,lambdaHat1 ))) ;

z = linspace(50,2000,200);
x = 0.5*(z(1:end-1)+z(2:end));
histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_power(x,lambdaHat1))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat1));

% MLE of power law with cut off
figure;
pf_powerC = @(x,alpha,lambda) x.^-alpha.*exp(-lambda*x)*lambda^(1-alpha)/gamma_incomplete(lambda*xmin,1-alpha) ;
[lambdaHat2,lambdaCI] = mle(sigIn, 'pdf',pf_powerC, 'start',[1,0], 'lowerbound',[1,0], 'upperbound', [3,10]) ;
L2 = sum(log(pf_powerC(sigIn,lambdaHat2(1), lambdaHat2(2) ))) ;


histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_powerC(x,lambdaHat2(1),lambdaHat2(2)))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_powerC(x, lambdaHat2(1),lambdaHat2(2)));


% MLE of power law with cut off
figure;
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*xmin) ;
[lambdaHat3,lambdaCI] = mle(sigIn, 'pdf',pf_exp, 'start',0, 'lowerbound',0, 'upperbound', 5) ;
L3 = sum(log(pf_exp(sigIn,lambdaHat3))) ;

histogram(sigIn,z,'Normalization','pdf');

hold on
plot(x,pf_exp(x,lambdaHat3))

figure
c = histcounts(sigIn,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));

%% generate distribution for exponential, power law and power law with cut off
%       x = randht(10000,'powerlaw',alpha);
%       x = randht(10000,'xmin',xmin,'powerlaw',alpha);
%       x = randht(10000,'cutoff',alpha, lambda);
%       x = randht(10000,'exponential',lambda);
%       x = randht(10000,'lognormal',mu,sigma);
%       x = randht(10000,'stretched',lambda,beta);

x = randht(length(patternScale), 'powerlaw', ) ;
