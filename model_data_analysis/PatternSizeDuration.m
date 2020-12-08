function c = PatternSizeDuration(mode)
close all;
clc;
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
hw = 31;
c = cell(2,num_files);
for j = 1:num_files
    Size = [];
    Duration = [];
    fprintf('Loading RYG.mat file %s...', files{j});
    R = load(files{j});
    R = get_grid_firing_centre(R);
    disp('done.\n');
    t_mid = R.grid.t_mid;
    ind_a_vec = R.grid.ind_ab(1,:);
    ind_b_vec = R.grid.ind_ab(2,:);
    
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
    
    [Lattice, ~] = lattice_nD(2, hw);
    x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
    y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
    
    i = 0;
    s = [];
    tic;
    for t = 1:length(t_mid)
        if ~isnan(x_centre(t))
            ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
            I = x_pos_o(ind_range_tmp);
            J = y_pos_o(ind_range_tmp);
            x_tmp = x_centre(t);
            y_tmp = y_centre(t);
            if i == 0
                x0 = x_tmp;
                y0 = y_tmp;
                i = i + 1;
                t0 = t_mid(t);
            end
            if Distance(x0,y0,x_tmp,y_tmp,2*hw + 1) > 18
                Size = [Size length(unique(s))];
                Duration = [Duration (t_mid(t)-t0)*0.1]; % ms
                s = [];
                t0 = t_mid(t);
            end
            ind = find(Distance(I,J,x_tmp,y_tmp,2*hw + 1) <= width(t));
            I = I(ind);
            J = J(ind);
            if length(I) < 1
                continue
            end
            s = [s IndexLattice(Lattice,I,J)];
            x0 = x_tmp;
            y0 = y_tmp;
        elseif ~isempty(s)
            Size = [Size length(unique(s))];
            Duration = [Duration (t_mid(t)-t0)*0.1]; % ms
            i = 0;
            s = [];
        end
    end
    c{1,j} = Size;
    c{2,j} = Duration;
    toc;
end
x = [];
y = [];
figure
for j = 1:num_files
    h = histogram(c{1,j},10);
    x = [x h.BinEdges(2:end)];
    y = [y h.Values/sum(h.Values)];
end
loglog(x,y,'.')
end