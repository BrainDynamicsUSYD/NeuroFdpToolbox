l = length(fre_dur);
Mat_distribution = zeros(l,400);
bins = [0:1:400];
for i = 1:l
    if length(fre_dur{i}) > 0
        fre_dur{i} = fre_dur{i}(fre_dur{i} > 3);
        Mat_distribution(i,:) = histcounts(fre_dur{i},bins);
    end
end
imagesc(Mat_distribution)
% imagesc(Mat_distribution(800:2401,:))
colorbar
set(gca,'XtickLabel',[50:50:400]); % ms
set(gca,'Ytick',[1:400:3601]);
set(gca,'YtickLabel',[100:(-10):10]);
% set(gca,'Ytick',[1:400:1601]);
% set(gca,'YtickLabel',[80:(-10):40]);
xlabel('duration (ms)')
ylabel('frequency(Hz)')