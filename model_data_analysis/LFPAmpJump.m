% jump behavior of LFP amplitude pattern
load('0006-201804201449-05982_in_1524200018506_0_neurosamp.mat','gamma_power_grid');
Y = prctile(gamma_power_grid(:),95);
s = size(gamma_power_grid);
mark = 0;
D = []; % propagate step size
L = []; % burst last time
xwidth = [];
ywidth = [];
finish = 1;
% subplot(1,3,1)
for t = 1:s(3)
    A = gamma_power_grid(:,:,t);
    [peak_mag,I] = max(A(:));
    [I_row, I_col] = ind2sub(s([1 2]),I);
    A_c = circshift(A, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
    GGrid = zeros(s([1 2]));
    GGrid(A_c >= Y) = 1;
    CC = bwconncomp(GGrid);
    validRegionIdx = 0;
    count = 1;
    for iRegion = 1:size(CC.PixelIdxList,2)
        if(size(CC.PixelIdxList{iRegion},1)>40)
            validRegionIdx(count) = iRegion ;
            count = count + 1 ;
        end
    end
    S = regionprops(CC,'Centroid');
    centroids = cat(1, S.Centroid);
    if count > 1
        centroids = centroids(validRegionIdx,:);
        centroids = centroids - [round(s(1)/2)-I_row  round(s(2)/2)-I_col];
        centroids = centroids - 32;
        if mark == 0
            temp = centroids;
            tempt = t;
        end
        if count == 2
            B = regionprops(CC,'BoundingBox'); % Smallest rectangle containing the region
            Boundary = cat(1, B.BoundingBox);
            xwidth = [xwidth Boundary(3)];
            ywidth = [ywidth Boundary(4)]; 
%             plot(centroids(1),centroids(2),'o')
%             hold on
            if t - tempt <= 20
                if finish >= 1
                    start = t;
                end
                finish = 0;
%                 plot([centroids(1),temp(1)],[centroids(2),temp(2)],'g')
                d = Distance_xy(centroids(1),centroids(2),temp(1),temp(2),s(1));
                D = [D d];
%                 hold on
            else
                finish = finish + 1;
                if finish == 1
                    l = tempt - start;
                    L = [L l];
                end
            end
            mark = mark + 1;
            temp = centroids;
            tempt = t;
        else
%             plot(centroids(:,1),centroids(:,2),'*')
%             hold on
        end
    end
    xlim([-31 31])
    ylim([-31 31])
end
D = 600/63*D; % um
L = 0.1*L; % ms
xwidth = 600/63*xwidth;  % um
ywidth = 600/63*ywidth;  % um
% text(-0.28,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(1,3,2)
histogram((xwidth+ywidth)/2)
xlabel('Width of Pattern(um)')
ylabel('Count')
text(-0.28,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
subplot(1,3,3)
histogram(D,1000)
xlim([0 20])
xlabel('Distance(um)')
ylabel('Count')
text(-0.28,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')