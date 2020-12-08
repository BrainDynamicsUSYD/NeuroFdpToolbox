figure_width = 11.4; % cm
figure_hight = 5.6; % cm
figure('NumberTitle','off','name', 'CH3Fig9', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

% R = load('0010-202001041836-25085_in_1578123572871_out_RYG2.mat');
% R = R.RL;
% load('3DBurstLFP0010minTime0SR1000P95.mat', 'distCent')
% load('3DBurstLFP0010minTime0SR1000P95.mat', 'WCentroids')

subplot(1,2,1)
histogram(R.grid.quick.jump_dist,60,'Normalization','probability')
xlabel('Jump Size(grid unit)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)

% distCent2 = [];
% for i = 1:length(WCentroids)
%     for j = 1:size(WCentroids{i},1)-1
%         distCent2 = [distCent2 Distance_xy(WCentroids{i}(j,2),WCentroids{i}(j,3),WCentroids{i}(j+1,2),WCentroids{i}(j+1,3),63)];
%     end
% end

subplot(1,2,2)
histogram(distCent,60,'Normalization','probability')
xlabel('Jump Size(grid unit)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.1,1,'B','Units', 'Normalized','FontSize',12)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig9 % this is the trick!!