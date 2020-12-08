%% This functions implements the joint distribution of peak amplitude and
% duration of the burst patterns as in the PNAS beta band paper
%
% Author: Xian Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% joint distribution
meanPower = sumAmp./patternScale ;
xVar = Duration' ;
yVar = peakAmp' ;
X = [xVar,yVar] ;
[N,C] = hist3(X,[50,1000]) ;
% xlabel('duration'); ylabel('mean power');
%%
Nlog= log(N./max(N(:))) ;
imagesc(C{1},C{2},Nlog)
set(gca,'YDir','normal')
% xlim([0,1600])

%%
% idx = find((Nlog(:)~=-Inf)) ;
% x = floor(idx/size(Nlog,1)) ;
% y = mod(idx,size(Nlog,1)) ;
% std = sqrt(1./Nlog(idx)*1000) ;
% data = [x, y, std] ;
% Results=weightedfit(data) ;

idx = find((N(:)~=0)) ;
x = floor(idx/size(N,1)) ;
y = mod(idx,size(N,1)) ;
std = sqrt(1./N(idx)) ;
data = [x, y, std] ;
Results=weightedfit(data) ;

%%
figure;
Nlog= log(N./max(N(:))) ;
imagesc(C{1},C{2},Nlog)
set(gca,'YDir','normal')
xlim([0 3300])
% xlim([0,1600])
cmap = parula ;
cmap(1,:) = ones(1,3) ;
colormap(cmap)
colorbar
title('Original data')
xlabel('Burst duration (s)')
ylabel('Peak power')

hold on
b = Results.slope;
a = Results.Intercept;
x = 1:size(N,2) ;
y_fit = (a + b*x-1)*(max(C{2})-min(C{2}))/size(N,1)+min(C{2});
plot(linspace(min(C{1}),max(C{1}),size(N,2)),y_fit,'k.--','linewidth',2)

%%
figure;
contourf(C{1},C{2},N)

%% joint distribution 2
meanPower = sumAmp./patternScale ;
xVar2 = Duration' ;
yVar2 = peakAmp' ;
X2 = [xVar2,yVar2] ;
[N2,C2] = hist3(X2,[1000,100]) ;
% xlabel('duration'); ylabel('mean power');
%%
figure;
Nlog2= log(N2./max(N2(:))) ;
imagesc(C2{1},C2{2},Nlog2)
set(gca,'YDir','normal')

%%
figure;
contourf(C2{1},C2{2},N2)

%%
% idx = find((Nlog2(:)~=-Inf)) ;
% x = floor(idx/size(Nlog,2)) ;
% y = mod(idx,size(Nlog,2)) ;
% std = Nlog2(idx).^2 ;
% data2 = [x, y, std] ;
% Results2=weightedfit(data2) ;

idx2 = find((N2(:)~=0)) ;
x2 = floor(idx2/size(N2,1)) ;
y2 = mod(idx2,size(N2,1)) ;
std2 = sqrt(1./N2(idx2)) ;
data2 = [x2, y2, std2] ;
Results2=weightedfit(data2) ;

%%
figure;
Nlog2= log(N2./max(N2(:))) ;
imagesc(C2{1},C2{2},Nlog2)
set(gca,'YDir','normal')
cmap = parula ;
cmap(1,:) = ones(1,3) ;
colormap(cmap)
colorbar
title('Original data')
xlabel('Burst duration (s)')
ylabel('Peak power')
% xlim([0,1600])
hold on
b2 = Results2.slope;
a2 = Results2.Intercept;
x2 = 1:size(N2,2) ;
y_fit2 = (a2 + b2*x2-1)*(max(C2{2})-min(C2{2}))/size(N2,1)+min(C2{2});
plot(linspace(min(C2{1}),max(C2{1}),size(N2,2)),y_fit2,'k.--','linewidth',2)