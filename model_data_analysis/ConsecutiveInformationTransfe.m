% function ConsecutiveInformationTransfer(varargin)
tic;
clear all;
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%
% r = load('0001-201701231616-50372_in_1485148826761_out_RYG.mat');
% % Loop number for PBS array job
% loop_num = 0;
% % center = round(rand*3969);
% center = 1535;
% fprintf(num2str(center));
% hw = 31;
% num = 10;
% moments = 200;
% [Lattice, ~] = lattice_nD(2, hw);
% Sam0 = FindNeurons(Lattice,hw,num,center);
% [dist] = lattice_nD_find_dist(Lattice,hw,center);
% [d,ind] = sort(dist);
% % choose the distance between the center and the other samples [23,23.03)
% Samcenter = ind(d>=23 & d<23.03);
% samnum = length(Samcenter);
% Sample = zeros(samnum,num);
% te = zeros(1,samnum);
% % style = ones(1,samnum);
% for i = 1:samnum
%     Sample(i,:) = FindNeurons(Lattice,hw,num,Samcenter(i));
% end
%
% for i = [-1:(-1):(-98)]
%     Opt.tauy = i;
%     loop_num = loop_num + 1;
%
%     % For PBS array job
%     if nargin ~= 0
%         PBS_ARRAYID = varargin{1};
%         if loop_num ~=  PBS_ARRAYID
%             continue;
%         end
%     end
%     Opt.taux = -1; % arbitrary value>= tauy
%     Opt.trperm = 4;
%     Opt.method = 'dr';
%     Opt.bias = 'qe';
%     toc;
%     fprintf('Calculating transfer entropy now...\n');
%     for i = 1:samnum
%         ntesh = [];
%         for j = 1e4:10:(r.step_tot - 1e3)
%             try
%                 [NTEsh] = transferentropy(r.spike_hist{1}(Sam0',j:(j+moments-1))',r.spike_hist{1}(Sample(i,:),j:(j+moments-1))',Opt,'NTEsh');
%                 NTEsh = NTEsh(NTEsh ~= Inf);
%                 NTEsh = NTEsh(NTEsh ~= -Inf);
%                 ntesh = [ntesh,nanmean(NTEsh)];
%             catch
%             end
%         end
%         te(i) = nanmean(ntesh);
%         toc;
%     end
%     save([num2str(center),'_',num2str(Opt.tauy),'transferE.mat'],'te','Opt')
%     fprintf('All finished.\n');
% %     x = [Lattice([center,Samcenter'],1)];
% %     y = [Lattice([center,Samcenter'],2)];
% %     plot(x(1),y(1),'r>',x(2:end),y(2:end),'bo');
% %     xlim([-31 31]);
% %     ylim([-31 31]);
% %     hold on;
% %     for i = 1:samnum
% %         if te(i) > 0
% %             style = '-';
% %         else
% %             style = '--';
% %         end
% %         plot1 = plot(x([1,i+1]),y([2,i+1]),style,'LineWidth',abs(10*te(i)/mean(te)));
% %         plot1.Color(4) = 0.5;
% %     end
% end
% end

load('0001-201707201244-31484_config_data.mat');
X = round(LFP_centre_x);
Y = round(LFP_centre_y);
hw = 31;
[Lattice, ~] = lattice_nD(2,hw);
v = 31.7; % um/ms = mm/s
Opt.taux = -1; % arbitrary value>= tauy
Opt.trperm = 4;
Opt.method = 'dr';
Opt.bias = 'qe';
% for i = 1:16
%     EN(i,:) = FindNeurons(Lattice,hw,10,X(i),Y(i));
% end
te = zeros(1,num_files);
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    R = get_grid_firing_centre(R);
    disp('done.\n');
    ntesh = 0;
    for repeat = 1:10
        j = 0;
        t_mid = R.grid.t_mid;
        ind_a_vec = R.grid.ind_ab(1,:);
        ind_b_vec = R.grid.ind_ab(2,:);
        x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
        y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);
        x_centre = R.grid.quick.centre(1,:);
        y_centre = R.grid.quick.centre(2,:);
        for t = 500:length(t_mid)
            if ~isnan(x_centre(t))
                ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
                x_tmp = x_centre(t);
                y_tmp = y_centre(t);
                rx = x_pos_o(ind_range_tmp);
                ry = y_pos_o(ind_range_tmp);
                ind = randperm(length(rx),10);
                for k = 1:10
                    receive(k) = find(Lattice(:,1) == rx(ind(k)) & Lattice(:,2) == ry(ind(k)));
                end
                try
                    Opt.tauy = -round(Distance(x_tmp,y_tmp,cx,cy,hw*2 + 1)*10/v*10); % 0.1ms
                    %             Opt.tauy = -(10 + j*10);
                    [NTEsh] = transferentropy(R.spike_hist{1}(send,(t_mid(t) - 19 - j*10):t_mid(t))',R.spike_hist{1}(receive,(t_mid(t) - 19 - j*10):t_mid(t))',Opt,'NTEsh');
                    NTEsh = NTEsh(NTEsh ~= Inf);
                    NTEsh = NTEsh(NTEsh ~= -Inf);
                    if ~isnan(nanmean(NTEsh))
                        ntesh = ntesh + nanmean(NTEsh);
                    end
                    %             fprintf('Calculating transfer entropy now...\n');
                catch
                end
                send = receive;
                cx = x_tmp;
                cy = y_tmp;
                j = 0;
            else
                j = j + 1;
            end
        end
        clear send cx cy
    end
    te(i) = ntesh/repeat;
    toc;
end
save('te.mat','te')
toc;