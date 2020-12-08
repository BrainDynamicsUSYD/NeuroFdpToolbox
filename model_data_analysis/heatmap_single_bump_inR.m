function heatmap_single_bump_inR(R,save_figure)
% weighted matrix reflects the inverse radius of the pattern(recorded every 1ms)
if nargin <= 1
    save_figure = 1; % default
end
if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end
hw = 31;
fw = 2*hw+1;
ind1 = find(~isnan(R.grid.quick.centre(1,:)));
ind2 = find(~isnan(R.grid.quick.centre(2,:)));
ind3 = find(~isnan(R.grid.quick.radius));
ind = intersect(ind1,ind2);
ind = intersect(ind,ind3);
x_mean_chosen = R.grid.quick.centre(1,ind);
y_mean_chosen = R.grid.quick.centre(2,ind);
dist_std_chosen = R.grid.quick.radius(ind);
Mat = zeros(fw+1);
time = datestr(now,'yyyymmddTHHMMSS');
if isempty(ind)
    imagesc(zeros(fw))
else
    for t = 1:length(ind)
        x_tmp = round(x_mean_chosen(t) + hw + 1);
        y_tmp = round(y_mean_chosen(t) + hw + 1);
        Mat(x_tmp,y_tmp) = Mat(x_tmp,y_tmp) + 1/dist_std_chosen(t);
    end
    Mat(1,1:fw) = Mat(1,1:fw) + Mat(fw+1,1:fw);
    Mat(1:fw,1) = Mat(1:fw,1) + Mat(1:fw,fw+1);
    Mat(1,1) = Mat(1,1) + Mat(fw+1,fw+1);
    Mat = Mat(1:fw,1:fw);
    Mat = permute(Mat,[2 1]);
    Mat = flipud(Mat);
    imagesc(Mat)
end
colorbar
sca = numel(find(Mat));
fprintf('loop:%04i scatter:%d\n',R.ExplVar.loop_num,sca);
tit = sprintf('loop number = %04i time:%s',R.ExplVar.loop_num,time);
name = sprintf('%04i_heatmap_%s.pdf',R.ExplVar.loop_num,time);
title(tit)
set(gcf,'renderer','zbuffer');
if save_figure == 1
    fprintf('\t Saving figure...');
    saveas(gca,name);
    fprintf('done.\n');
else
    fprintf('\t Figure is not saved!\n');
end
end