% function AvalancheOnSpikes % (R)
dir_strut = dir('*_out_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% for i = 1:num_files
%     fprintf('Loading RYG.mat file %s...\n', files{i});
%     R = load(files{i});
% end

% %% mean ISI as the bin
% vec = sum(R.spike_hist{1});
% ind = find(vec);
% ISI = ind(2:end)-ind(1:end-1);
% mean(R.dt*ISI) % ms
%%

j = 1;
for bin = 0.1 % [0.2 1 2 4] % [1 2 4 8 16] % ms
    Size = [];
    Duration = [];
%     for n = [7:13:num_files] % 1:num_files
%         fprintf('Loading RYG.mat file %s...\n', files{n});
%         R = load(files{n});
        vec = full(sum(R.spike_hist{1}));
        %     mat = mat(1:end-mod(length(mat),bin/R.dt));
        mat = sum(reshape(vec,bin/R.dt,[]),1);
        NonA = find(mat==0);
        for i = 1:length(NonA)-1
            if NonA(i+1) - NonA(i) == 1
                continue
            end
            Size = [Size sum(mat(NonA(i):NonA(i+1)))];
            Duration = [Duration (NonA(i+1)-NonA(i)-1)*bin]; % ms
        end
%     end
    %     Scell{j} = Size;
    %     Dcell{j} = Duration;
    %     mean(Duration)
    %     Duration = Duration/bin;
    %     %% Plotting & Fitting
    %     bin = [0.2 1 2 4];
    %     for j = 1:4
    %         Size = Scell{j};
    %         Duration = Dcell{j};
    subplot(1,2,1) % (4,2,2*j-1)
    %     [N,edges] = histcounts(Duration,200);
    edges = 10.^(linspace(log10(min(Size)),log10(max(Size)),6001));
    N = histcounts(Size,edges);
    Y = (edges(1:end-1)+edges(2:end))/2;
    loglog(Y,N,'o')
    hold on
    Y = Y(N > 0);
    N = N(N > 0);
    %     Y = Y(1:end-round(end/2));
    %     N = N(1:end-round(end/2));
    v = polyfit(log10(Y),log10(N),1);
    x = Y;
    y = 10^v(2)*x.^v(1);
    loglog(x,y)
    str = ['p = ',num2str(v(1))];
    text(min(x),max(y),str)
    xlabel('Size')
    ylabel('Count')
    legend(strcat('bin=',num2str(bin),'ms'))
    subplot(1,2,2) % (4,2,j*2)
    %     [N,edges] = histcounts(Duration,60);
    edges = 10.^(linspace(log10(min(Duration)),log10(max(Duration)),6001));
    N = histcounts(Size,edges);
    Y = (edges(1:end-1)+edges(2:end))/2;
    loglog(Y,N,'o')
    hold on
    Y = Y(N > 0);
    N = N(N > 0);
    %     Y = Y(1:end-round(end/2));
    %     N = N(1:end-round(end/2));
    v = polyfit(log10(Y),log10(N),1);
    x = Y;
    y = 10^v(2)*x.^v(1);
    loglog(x,y)
    str = ['p = ',num2str(v(1))];
    text(min(x),max(y),str)
    xlabel('Duration')
    ylabel('Count')
    legend(strcat('bin=',num2str(bin),'ms'))
end
% %     %%
%     j = j + 1;

%% TOOLBOX
% pdf: y = x^-alpha
% cdf: slope is -alpha+1
[alpha, xmin, L]=plfit(Size)
%             xmin = 4; alpha = 2.88;
plplot(Size, xmin, alpha);
% [p,gof]=plpva(patternScale, xmin)
str = ['Size pl:alpha = ',num2str(alpha)];
title(str)

[alpha, xmin, L]=plfit(Duration)
% xmin = 0.2; alpha = 3.7;
plplot(Duration,xmin,alpha);
% [p,gof]=plpva(patternScale, xmin)
str = ['Duration pl:alpha = ',num2str(alpha)]; % ,'(MF)'];
title(str)

%% XIAN'S COMPRISON BETWEEN POWER LAW AND EXPONENTIAL
pf_power = @(x,alpha,a) x.^-alpha*(alpha-1)*a^(alpha-1);
[lambdaHat,lambdaCI] = mle(Size, 'pdf',pf_power, 'start',[1.1,3], 'lowerbound',[1,1], 'upperbound', [5,745]) ;
L = sum(log(pf_power(Size,lambdaHat(1), lambdaHat(2) )))
lambdaHat
pf_exp = @(x,lambda,a) exp(-lambda*x) *lambda*exp(lambda*a) ;
[lambdaHat3,lambdaCI] = mle(Size, 'pdf',pf_exp, 'start',[0,3] , 'lowerbound',[0,1], 'upperbound', [5,745]) ;
L3 = sum(log(pf_exp(Size,lambdaHat3(1),lambdaHat3(2))))
lambdaHat3

pf_power = @(x,alpha,a) x.^-alpha*(alpha-1)*a^(alpha-1);
[lambdaHat,lambdaCI] = mle(Duration, 'pdf',pf_power, 'start',[1.1,0.2], 'lowerbound',[1,0.1], 'upperbound', [5,34]) ;
L = sum(log(pf_power(Duration,lambdaHat(1), lambdaHat(2) )))
lambdaHat
pf_exp = @(x,lambda,a) exp(-lambda*x) *lambda*exp(lambda*a) ;
[lambdaHat3,lambdaCI] = mle(Duration, 'pdf',pf_exp, 'start',[0,0.2] , 'lowerbound',[0,0.1], 'upperbound', [5,34]) ;
L3 = sum(log(pf_exp(Duration,lambdaHat3(1),lambdaHat3(2))))
lambdaHat3

%% DEFINE a(starting value) already
signal = Size; % (Size>2890);
a = 3 ;
pf_power = @(x,alpha) x.^-alpha*(alpha-1)*a^(alpha-1);
[lambdaHat,lambdaCI] = mle(signal, 'pdf',pf_power, 'start',[1.1], 'lowerbound',[1], 'upperbound', [3.5]) ;
L = sum(log(pf_power(signal,lambdaHat(1))))
lambdaHat

% % exponential distribution
pf_exp = @(x,lambda) exp(-lambda*x) *lambda*exp(lambda*a) ;
[lambdaHat3,lambdaCI] = mle(signal, 'pdf',pf_exp, 'start',[0] , 'lowerbound',[0], 'upperbound', [2]) ;
L3 = sum(log(pf_exp(signal,lambdaHat3(1))))
lambdaHat3
%%
z = linspace(0.2,34,2e2); % patternScale:580,9210; Duration:30,100,50
x = 0.5*(z(1:end-1)+z(2:end));
subplot(2,2,1)
histogram(signal,z,'Normalization','pdf')
hold on
plot(x,pf_power(x,lambdaHat))
subplot(2,2,2)
c = histcounts(signal,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_power(x, lambdaHat));
str = ['p = ',num2str(-lambdaHat),',L = ',num2str(L)];
title(str)
subplot(2,2,3)
histogram(signal,z,'Normalization','pdf')
hold on
plot(x,pf_exp(x,lambdaHat3))
xlabel('Size')
subplot(2,2,4)
c = histcounts(signal,z,'Normalization','pdf');
loglog(x,c,'.','markersize',20)
hold on
loglog(x,pf_exp(x, lambdaHat3));
xlabel('Size')
str = ['L = ',num2str(L3)];
title(str)

%% NCC TOOLBOX
[tau, xmin, xmax, L] = plmle(Size,'xmin',5); % ,'xmax',38);
fitParams = struct; fitParams.tau = tau;
fitParams.xmin = xmin; % xmin; % min(Size);
fitParams.xmax = xmax; % xmax; % max(Size);
plotParams = struct; plotParams.dot = 'on';
titleString = ['Size pl:alpha = ',num2str(tau,'%.4f')];
data = plplottool(Size,'plotParams',plotParams,'fitParams',fitParams, 'title',titleString);
%%
[tau, xmin, xmax, L] = plmle(Duration,'xmin',0.3,'xmax',1.2); % ,'xmin',0.3,'xmax',1.2); % ) % 
fitParams = struct; fitParams.tau = tau;
fitParams.xmax = max(Duration); % max(Duration); % xmax; % max(Duration);
fitParams.xmin = xmin; % min(Duration); % xmin; % min(Duration);
plotParams = struct; plotParams.dot = 'on';
titleString = ['Duration pl:alpha = ',num2str(tau,'%.4f')];
data = plplottool(Duration,'plotParams',plotParams,'fitParams',fitParams, 'title',titleString);
% end
% end