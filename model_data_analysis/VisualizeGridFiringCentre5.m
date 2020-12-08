function VisualizeGridFiringCentre5(stimulus,ind)
% adapt from function VisualizeGridFiringCentre2.m
% only show case which meet my requirement
% optional: adding local hotpot
% close all;
clc;
mode = 'quick';
hotpot = [513 2497];
hw = 31;
fw = 2*hw+1;
[Lattice, ~] = lattice_nD(2, hw);
post_dist = lattice_nD_find_dist(Lattice,hw,1);
[post_dist,~] = sort(post_dist);
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
stimulus(1) = stimulus(1) - 20; % leave some time space
stimulus(end) = stimulus(end) + 20; % leave some time space
otherind = true(1,length(hotpot));
otherind(ind) = false;
for id_out = 151:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    R = get_grid_firing_centre(R);
    t_mid = R.grid.t_mid;
    [~,p1] = sort(abs(t_mid-stimulus(1)));
    p1 = p1(1);
    [~,p2] = sort(abs(t_mid-stimulus(end)));
    p2 = p2(1);
    r = post_dist(R.ExplVar.local_population);
    nump1 = zeros(1,length(hotpot));
    nump2 = nump1;
    switch mode
        case 'bayesian'
            x_centre = R.grid.bayes.centre(1,:);
            y_centre = R.grid.bayes.centre(2,:);
            width = R.grid.bayes.radius;
        case 'quick'
            x_centre = R.grid.quick.centre(1,:);
            y_centre = R.grid.quick.centre(2,:);
            width = R.grid.quick.radius;
    end
    for i = 1:length(hotpot)
        d = Distance_xy(x_centre(1:p1),y_centre(1:p1),Lattice(hotpot(i),1),Lattice(hotpot(i),2),fw);
        nump1(i) = sum(d <= r);
        d = Distance_xy(x_centre(p2:end),y_centre(p2:end),Lattice(hotpot(i),1),Lattice(hotpot(i),2),fw);
        nump2(i) = sum(d <= r);
    end
    if max(nump1)-min(nump1)<max([0.5*max(nump1) 1]) && nump2(ind)>1.5*max(nump2(otherind)) && nump2(ind)/(pi*r^2) > 1.1*(length(x_centre)-p2)/fw^2
        subplot(1,3,1)
        for i = 1:length(hotpot)
            circle(Lattice(hotpot(i),1),Lattice(hotpot(i),2),r)
            hold on;
        end
        plot(x_centre(1:p1),y_centre(1:p1),'r.', 'MarkerSize', 8)
        ts = sprintf('time = %8.1f to %8.1fms', t_mid(1)*0.1,t_mid(p1)*0.1);
        title(ts);
        xlim([-hw hw]);
        ylim([-hw hw]);
        subplot(1,3,2)
        for i = 1:length(hotpot)
            circle(Lattice(hotpot(i),1),Lattice(hotpot(i),2),r)
            hold on;
        end
        plot(x_centre(p1:p2),y_centre(p1:p2),'r.', 'MarkerSize', 8)
        ts = sprintf('time = %8.1f to %8.1fms', t_mid(p1)*0.1,t_mid(p2)*0.1);
        title(ts);
        xlim([-hw hw]);
        ylim([-hw hw]);
        subplot(1,3,3)
        for i = 1:length(hotpot)
            circle(Lattice(hotpot(i),1),Lattice(hotpot(i),2),r)
            hold on;
        end
        plot(x_centre(p2:end),y_centre(p2:end),'r.', 'MarkerSize', 8)
        ts = sprintf('time = %8.1f to %8.1fms', t_mid(p2)*0.1,t_mid(end)*0.1);
        title(ts);
        xlim([-hw hw]);
        ylim([-hw hw]);
        next = input('\t Next figure?');
        close all
    end
end
end