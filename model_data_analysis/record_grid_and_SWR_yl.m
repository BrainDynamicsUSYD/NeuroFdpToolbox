function record_grid_and_SWR_yl(R, LFP_centre_x, LFP_centre_y)
close all;
clc;
hw = 31;
dt = R.dt;
t_mid = R.grid.raw.t_mid;
figure('Name','Vis','color','w','NumberTitle','off');
axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);
xlim([-hw hw]);
ylim([-hw hw]);
text(LFP_centre_x(:),LFP_centre_y(:),{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'; ...
    '11';'12';'13';'14';'15';'16'});
%     '11';'12';'13';'14';'15';'16';'17';'18';'19';'20'; ...
%     '21';'22';'23';'24';'25';'26';'27';'28';'29';'30'; ...
%     '31';'32';'33';'34';'35';'36';'37';'38';'39';'40'; ...
%     '41';'42';'43';'44';'45';'46';'47';'48';'49';'50'; ...
%     '51';'52';'53';'54';'55';'56';'57';'58';'59';'60'; ...
%     '61';'62';'63';'64'});
hold on;
fid = fopen('SWR_timepoint0001B.txt','w');
for tt = unique(sort([1:10:t_mid(end) t_mid(:)']))
    h2 = plot(100,0);
    h3 = plot(100,0);
    % ripple dectection
    detceded = sum(R.LFP.ripple_event.is_SWR(:, tt), 2) > 0;
    ind = find(detceded);
    if length(ind) == 0
        h4 = plot(LFP_centre_x(detceded), LFP_centre_y(detceded), '.r', 'MarkerSize',6);
    else
        mat_d = GridDM([4 4 1], [4 4 1], 1);
        prin_d = mat_d(ind(:),ind(:)) > sqrt(2);
        prin_group = sum(prin_d,2) == (length(ind) - 1);
        ind_group = find(prin_group);
        ind_combine = find(~prin_group);
        prin_groups = sum(prin_group);
        switch prin_groups
            case 4
                h4 = gscatter(LFP_centre_x(detceded), LFP_centre_y(detceded), [1:prin_groups]);
                fprintf(fid,[ts]);
                fprintf(fid,'(%d)\n',prin_groups);
            case 3
                h4 = gscatter(LFP_centre_x(detceded), LFP_centre_y(detceded), [1:prin_groups]);
                fprintf(fid,[ts]);
                fprintf(fid,'(%d)\n',prin_groups);
            case 2
                if length(ind) == 2
                    h4 = gscatter(LFP_centre_x(detceded), LFP_centre_y(detceded), [1:prin_groups]);
                    fprintf(fid,[ts]);
                    fprintf(fid,'(%d)\n',prin_groups);
                else
                    x_combine = LFP_centre_x(ind(ind_combine));
                    y_combine = LFP_centre_y(ind(ind_combine));
                    [x,y] = combination(x_combine,y_combine,hw);
                    x = [x;LFP_centre_x(ind(ind_group))];
                    y = [y;LFP_centre_y(ind(ind_group))];
                    h4 = gscatter(x, y, [1:prin_groups + 1]);
                    fprintf(fid,[ts]);
                    fprintf(fid,'(%d)\n',prin_groups + 1);
                end
            case 1
                x_combine = LFP_centre_x(ind(ind_combine));
                y_combine = LFP_centre_y(ind(ind_combine));
                [x,y] = combination(x_combine,y_combine,hw);
                x = [x;LFP_centre_x(ind(ind_group))];
                y = [y;LFP_centre_y(ind(ind_group))];
                h4 = gscatter(x, y, [1:length(x)]);
                fprintf(fid,[ts]);
                fprintf(fid,'(%d)\n',length(x));
            case 0
                ind_first1 = find(prin_d(1,:));
                ind_first0 = find(~prin_d(1,:));
                ind_1group = ind(ind_first1(logical(prod(prin_d(ind_first1,ind_first0),2))));
                x_combine1 = LFP_centre_x(ind_1group);
                y_combine1 = LFP_centre_y(ind_1group);
                x_combine2 = LFP_centre_x(ind(~ismember(ind,ind_1group)));
                y_combine2 = LFP_centre_y(ind(~ismember(ind,ind_1group)));
                [x1,y1] = combination(x_combine1,y_combine1,hw);
                [x2,y2] = combination(x_combine2,y_combine2,hw);
                x = [x1;x2];
                y = [y1;y2];
                h4 = gscatter(x, y, [1:length(x)]);
                fprintf(fid,[ts]);
                fprintf(fid,'(%d)\n',length(x));
        end
        fprintf(fid,'Detector:%d ',ind);
        fprintf(fid,'\n');
    end
    pause(0.1);
    delete(h2);
    delete(h3);
    delete(h4);
    ts = sprintf('time = %8.1f ms', tt*dt);
    title(ts);
end
fclose(fid);
end
