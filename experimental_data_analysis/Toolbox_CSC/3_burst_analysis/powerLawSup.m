%% scale vs duration
[sortDuration,durataionIdx] = sort(Duration) ;
sortScale = patternScale(durataionIdx) ;
[DuSort,~,DuIdx] = unique(Duration) ;
ScaleDu = zeros(1,length(DuSort)) ;
idxCell = cell(length(DuSort),1) ;

for k = 1:length(DuSort)
    idxCell{k} = find(DuIdx == k) ;
    ScaleDu(k) = mean(patternScale(idxCell{k})) ;
end


loglog(Duration,patternScale,'o') ;
figure;
loglog(DuSort,ScaleDu,'.','markersize',12) 

p = polyfit(log(DuSort),log(ScaleDu),1)

alpha = plfit(Duration) ;
tau = plfit(patternScale) ;

(alpha - 1)/(tau -1 )
p

%% scale profile
for iBurst = 1:size(instantScale,2)
    plot(instantScale{iBurst})
    xlim([0,1000])
    ylim([0 50])
    %pause
    %close all
    hold on
    
end


%% scale profile
[sortDuration,durataionIdx] = sort(Duration) ;
sortScale = patternScale(durataionIdx) ;
[DuSort,~,DuIdx] = unique(Duration) ;
idxCell = cell(length(DuSort),1) ;

meanInstantScale = zeros(4000, length(DuSort)) ;
for k = 1:length(DuSort)
    idxCell{k} = find(DuIdx == k) ;
    tempInstantScale = zeros(4000,1) ;
    for j = 1:length(idxCell{k})
        temp = idxCell{k} ;
        tempInstantScale = tempInstantScale+[instantScale{temp(j)};...
            zeros(length(tempInstantScale)-length(instantScale{temp(j)}),1) ];
    end
    meanInstantScale(:,k) = tempInstantScale/j ;

end
% start of indice
 DuSort(1)
 % end of indice
  DuSort(end)
%%
figure;
plot(meanInstantScale(:,6))
hold on
plot(meanInstantScale(:,15))
hold on
plot(meanInstantScale(:,25))
hold on
plot(meanInstantScale(:,35))
hold on
plot(meanInstantScale(:,45))
hold on
plot(meanInstantScale(:,51))
legend('30ms Duration','40ms Duration', '50ms Duration','60ms Duration',...
    '70ms Duration','80ms Duration')

%%
figure;
x = 1:4000 ;
plot(x/30*80,meanInstantScale(:,25)*30^(1-p(1)))
hold on
plot(x/40*80,meanInstantScale(:,35)*40^(1-p(1)))
hold on
plot(x/50*80,meanInstantScale(:,45)*50^(1-p(1)))
hold on
plot(x/60*80,meanInstantScale(:,55)*60^(1-p(1)))
hold on
plot(x/70*80,meanInstantScale(:,65)*70^(1-p(1)))
hold on
plot(x/80*80,meanInstantScale(:,75)*80^(1-p(1)))
hold on
% plot(x,meanInstantScale(:,51)*80^(1-p(1)))
legend('30ms Duration','40ms Duration', '50ms Duration','60ms Duration',...
    '70ms Duration','80ms Duration')