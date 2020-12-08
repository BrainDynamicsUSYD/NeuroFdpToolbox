is_SWR = R.LFP.ripple_event.is_SWR(:,R.grid.t_mid);
ind = find(sum(is_SWR) > 0);
is_SWR = logical(is_SWR(:,ind));
centre = R.grid.centre(:,ind);
radius = R.grid.radius(ind);
t_match = R.grid.t_mid(ind);
fr_E = R.Analysis.Hz_t{1};
fr_E = fr_E(t_match);
n_content = length(ind);
l = 62;
d = zeros(1,n_content);
for t = 1:n_content
    I = LFP_centre_x(is_SWR(:,t));
    J = LFP_centre_y(is_SWR(:,t));
    x = centre(1,t);
    y = centre(2,t);
    d(t) = mean(Distance(I,J,x,y,l));
end
d_r = d./radius;
figure(1)
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Distance between Bump and SWRs', 'FontSize', 12, 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
subplot(1,2,1)
histogram(d,25)
xlabel({['Distance,mean=',num2str(mean(d))],['max=',num2str(max(d)),',min=',num2str(min(d))]},'FontSize',10)
subplot(1,2,2)
histogram(d_r,25)
xlabel({['Normalized Distance,mean=',num2str(mean(d_r))],['max=',num2str(max(d_r)),',min=',num2str(min(d_r))]},'FontSize',10)
figure(2)
histogram(fr_E,25)
title('Distribution of Bump Firing Rate')
xlabel(['mean=',num2str(mean(fr_E)),'Hz,min=',num2str(min(fr_E)),'Hz,max=',num2str(max(fr_E)),'Hz'])