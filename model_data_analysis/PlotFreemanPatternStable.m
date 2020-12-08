function q = PlotFreemanPatternStable(ApiceP,ApiceN)
% example on FreemanPatternStable Aone/two
% Ref: Freeman, Walter J., and John M. Barrie. "Analysis of spatial patterns
%     of phase in neocortical gamma EEGs in rabbit." Journal of neurophysiology 84.3 (2000): 1266-1278.

%% calculate q vector
qP = zeros(1,length(ApiceP)); % positive cone q vector
i = 1;
while i <= length(ApiceP)
    qP(i) = size(ApiceP{i},2)/2 -2;
    if qP(i) < 0
        i = i + 1;
    else
        qP(i+1:i+qP(i)+1) = 0;
        i = i + qP(i) + 2;
    end
end
qN = zeros(1,length(ApiceN)); % negative cone q vector
i = 1;
while i <= length(ApiceN)
    qN(i) = size(ApiceN{i},2)/2 -2;
    if qN(i) < 0
        i = i + 1;
    else
        qN(i+1:i+qN(i)+1) = 0;
        i = i + qN(i) + 2;
    end
end
q = qP;
Apice = ApiceP;
IndNq = find(qN >= 25);
IndPminus = find(q < 0);
cover1 = IndPminus(ismember(IndPminus,IndNq));
CoverInd1 = cover1(cover1 > 0); % index of q>= 25 in Nvector covered in Pvector
q(CoverInd1) = qN(CoverInd1);
Apice(CoverInd1) = ApiceN(CoverInd1);
for i = 1:length(CoverInd1)
    ind = CoverInd1(i);
    StopInd = find(q(ind+1:ind+q(ind)+1)>0);
    if isempty(StopInd)
        q(ind+1:ind+q(ind)+1) = 0;
        Apice(ind+1:ind+q(ind)+1) = ApiceN(ind+1:ind+q(ind)+1);
    else
        q(ind+1:ind+StopInd(1)-2) = 0;
        Apice{ind} = Apice{ind}(:,1:2*StopInd(1));
        Apice(ind+1:ind+StopInd(1)-2) = ApiceN(ind+1:ind+StopInd(1)-2);
    end
end
IndNminusOne = find(qN == -1);
IndPminusTwo = find(q == -2);
cover2 = IndPminusTwo(ismember(IndPminusTwo,IndNminusOne));
CoverInd2 = cover2(cover2 > 0);
q(CoverInd2) = -1;
Apice(CoverInd2) = ApiceN(CoverInd2);

subplot(3,2,6)
%% plot stable segments figure
line([0 10 10],[0 0 10],'Color','black','LineStyle','-')
hold on;
line([0 0 10],[0 10 10],'Color','black','LineStyle','-')
xlim([-10 20])
ylim([-10 20])
l = length(Apice);
i = 1990;
pre = [];
while i <= 2068 % l
    if q(i) == -2
        delete(findobj(gca,'Type','line','Color','g'));
        delete(findobj(gca,'Color','r'));
        delete(findobj(gca,'Color','b'));
        delete(findobj(gca,'Color','y'));
        pre = [];
        i = i + 1;
    else
        if isempty(pre)
            pre = Apice{i}(1,1:2);
        end
        if q(i) == -1
            vec = Apice{i};
            hold on
            plot(vec(1),vec(2), 'b.', 'MarkerSize', 4);
            hold on
            plot([pre(1),vec(1)],[pre(2),vec(2)],'y')
            pre = vec;
            i = i + 1;
        else
            vec = Apice{i}(1,:);
            plot(vec(1:2:end),vec(2:2:end), 'r.', 'MarkerSize', 8);
            hold on
            plot([pre(1),vec(1)],[pre(2),vec(2)],'g')
            hold on
            plot(vec(1:2:end),vec(2:2:end),'g')
            pre = vec;
            i = i + q(i) + 2;
        end
    end
    %     ts = sprintf('Time segment from %8.1f ms to %8.1f ms', (i-1)*2,(i-1)*2 +128);
    %     title(ts);
    text(-0.28,1.02,'F','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    pause(0.2)
end
end