function similarity_heatmap_connections(varargin)
d_RYG = dir('*RYG.mat');
d_ygin_syn= dir('*ygin_syn');
r = 7;
type1 = 'EE';
type2 = 'EIE';
h = 0;

% Loop number for PBS array job
loop_num = 0;

for m = 1:150
    loop_num = loop_num + 1;
    
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    tic;
    R = load(d_RYG(m).name);
    hw = R.grid.hw;
    fw = 2*hw+1;
    Mat = zeros(fw);
    time = datestr(now,'yyyymmddTHHMMSS');
    
    if size(R.grid.centre) == [0,0]
        imagesc(Mat)
    else
        x_mean_chosen = R.grid.centre(1,:);
        y_mean_chosen = R.grid.centre(2,:);
        dist_std_chosen = R.grid.radius;
        t_mid = R.grid.t_mid_full;
        t_mid_chosen = R.grid.t_mid;
        j = 1;
        for t = 1:length(t_mid)
            if sum(t_mid_chosen == t_mid(t)) == 1
                x_tmp = round(x_mean_chosen(j)) + hw + 1;
                y_tmp = round(y_mean_chosen(j)) + hw + 1;
                Mat(x_tmp,y_tmp) = Mat(x_tmp,y_tmp) + 1/dist_std_chosen(j);
                j = j + 1;
            end
        end
        Mat = permute(Mat,[2 1]);
        Mat = flipud(Mat);
        Mat = (Mat - mean(Mat(:)))/std(Mat(:));%standard normalized the matrix
        imagesc(Mat)
    end
    name_heatmap = sprintf('%04i_heatmap_%s.jpg',R.ExplVar.loop_num,time);
    set(gcf,'renderer','zbuffer');
    saveas(gcf,name_heatmap);
    AA = imread(name_heatmap);
    delete(name_heatmap);
    fprintf('Reading bump heatmap done...\n');
    
    %     Mat_fr = flipud(vec2mat(R.Analysis.rate{1},fw,fw));
    %     Mat_fr = (Mat_fr - min(Mat_fr(:)))/(max(Mat_fr(:)) - min(Mat_fr(:)));%normalized the matrix
    %     imagesc(Mat_fr);
    %     %     colorbar
    %     %     time = datestr(now,'yyyymmddTHHMMSS');
    %     %     tit = sprintf('firing rate image of excitatory neurons,loop number = %04i time:%s',R.ExplVar.loop_num,time);
    %     name_fr = sprintf('%04i_image_firing_rate_%s.jpg',R.ExplVar.loop_num,time);
    %     %     title(tit)
    %     set(gcf,'renderer','zbuffer');
    %     saveas(gcf,name_fr);
    %     BB = imread(name_fr);
    %     delete(name_fr);
    %     fprintf('Reading firing rate image done...\n');
    
    mat_EE = zeros(fw);
    mat_EIE = zeros(fw);
    [Lattice, N] = lattice_nD(2, hw);
    R_EE = read_ygin_syn(d_ygin_syn(m).name,1,1);
    R_IE = read_ygin_syn(d_ygin_syn(m).name,1,2);
    R_EI = read_ygin_syn(d_ygin_syn(m).name,2,1);
    R_EE = R_EE{1};
    R_IE = R_IE{1};
    R_EI = R_EI{1};
    for i = 1:N
        dist = lattice_nD_find_dist(Lattice, hw, i);
        local = (find(dist <= r));
        con_EE = ismember(R_EE.I,local).*ismember(R_EE.J,local);% EE
        IE = R_IE.J(ismember(R_IE.I,local));%find the connection from local E to the inhibitory population
        EI = R_EI.I(ismember(R_EI.J,local));%find the connection from the inhibitory population to local E
        EIE = intersect(IE,EI);%find the indices of inhibitory neurons which both receive connections from local E and send connection to local E
        
        %count the connections of individual EIE neuron receiving from local E
        [count1,element1] = hist(IE,unique(IE));
        num_pre = count1(ismember(element1,EIE));
        
        %count the connections of individual EIE neuron sending to local E
        [count2,element2] = hist(EI,unique(EI));
        num_post = count2(ismember(element2,EIE));
        
        num = sum(num_pre.*num_post);%get the total number of EIE-loop
        mat_EE(i) = sum(R_EE.K(logical(con_EE)));
        mat_EIE(i) = mean(R_IE.K)*mean(R_EI.K)*num;% the whole connection stengths of EIE-loop
    end
    mat_EE = (mat_EE - mean(mat_EE(:)))/std(mat_EE(:));
    mat_EIE = (mat_EIE - mean(mat_EIE(:)))/std(mat_EIE(:));
    mat_EIE = -1*mat_EIE;
    
    % the reason why we don't consider II effect is the possibility IE and
    % EI are relatively high.When the radius of local field is not less
    % than 3,the inhibitory neurons participate in the EIIE loop are the
    % whole Inhibitory population.So the effect would be a constant,which
    % make no difference on the local field connection images.
    
    fprintf('Work before chosing the coefficients done...\n');
    toc;
    
    for k = [0:0.01:1]
        mat_c_k = (1 - k)*mat_EE + k*mat_EIE;
        mat_c_k = permute(mat_c_k,[2 1]);
        mat_c_k = flipud(mat_c_k);
        mat_c_k = (mat_c_k - mean(mat_c_k(:)))/std(mat_c_k(:));%normalized the matrix
        figure(1)
        imagesc(mat_c_k)
        %colorbar
        set(gcf,'renderer','zbuffer');
        %         t2 = sprintf(['local field connections %s-%0.5g%s+%0.5g%s,loop number = %04i,r = %d'],type1,k1,type2,k2,type3,R.ExplVar.loop_num,r);
        %         title(t2);
        name_connections = sprintf(['%04i_local_field_connections_%0.2g%s-%0.2g%s_r%d.jpg'],m,1 - k,type1,k,type2,r);
        saveas(gca,name_connections);
        CC = imread(name_connections);
        delete(name_connections);
        if ssim(CC,AA) > h
            h = ssim(CC,AA);
            k_best_bump = k;
        end
    end
    fprintf('Best weight for firing rate:%0.2g\nSimilarity value:%0.4g\n',k_best_bump,h)
    toc;
end
end
