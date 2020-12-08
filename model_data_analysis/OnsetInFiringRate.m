% metastable point in firing rate
[fr,~] = CollectCellYG('Analysis','mean(Analysis.rate{1})');
fr = vec2mat([fr{:}],13);
Fr = mean(fr);
err = std(fr);
ratio = 0.7:0.05:1.3;
x1 = 0.7:0.01:1.05;
x2 = 0.9:0.01:1.3;
for i = 1:2
    figure;
    errorbar(ratio(1:7),Fr(1:7),err(1:7),'ro');
    hold on;
    errorbar(ratio(7:13),Fr(7:13),err(7:13),'b*');
    hold on;
    %     y1 = -278.92*x1.^3 + 1218*x1.^2 - 1602*x1 + 665.22;
    %     y2 = -49.264*x2.^3 + 181.75*x2.^2 - 225.42*x2 + 95.261;    
    if i == 2
        y1 = exp(-7.1088*x1.^2 + 1.7663*x1 + 6.0946);
        y2 = exp(4.5*x2.^2 - 12.643*x2 + 8.9733);
        set(gca, 'YScale', 'log');
        plot(x1,y1,'r',x2,y2,'b--')
    else
        y1 = exp(-7.1088*ratio.^2 + 1.7663*ratio + 6.0946);
        y2 = exp(4.5*ratio.^2 - 12.643*ratio + 8.9733);
        plot(ratio,y1,'r',ratio,y2,'b--')
    end
    xlabel('F_r')
    ylabel('Firing rate(Hz)')
    % set(gcf,'color','w');
end