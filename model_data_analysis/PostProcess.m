function PostProcess % (varargin)
% run identical analysis on all RYG.mat on cluster
%%
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
dir_strut = dir('*_RYG.mat');
xlsfiles={dir_strut.name};
[~,idx]=natsortfiles(xlsfiles);
dir_strut=dir_strut(idx);
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
% dir_strut2 = dir('*_config_data.mat');
% num_files2 = length(dir_strut2);
% files2 = cell(1,num_files2);
% for id_out = 1:num_files2
%     files2{id_out} = dir_strut2(id_out).name;
% end
% load('StiNeu.mat')
% fs = 1e3;
% dt = 1;
% lowFreq = 30;
% hiFreq = 80; % Hz
% Wn = [lowFreq hiFreq]/(fs/2);
% order = 4; % 4th order
% [b,a] = butter(order/2,Wn,'bandpass');
% gaus_width = 12.5; % ms
% [ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );

% H = [];
% L = [];
% series = zeros(1,13);
% start = [8.5e4 10.5e4 12.5e4 5.5e4 14.5e4 2.5e4*ones(1,2) 5.5e4*ones(1,2) 2.5e4*ones(1,2) 6.5e4 2.5e4];
% i = 1;
% % Loop number for PBS array job
% loop_num = 0;
% for id_out = 1:num_files
%     % For PBS array job
%     loop_num = loop_num + 1;
%     if nargin ~= 0
%         PBS_ARRAYID = varargin{1};
%         if loop_num ~=  PBS_ARRAYID
%             continue;
%         end
%     end
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
%     %     power = mean(R.LFP.LFP_gamma_hilbert_abs);
%     %     series(i) = mean(power(start(i):start(i)+4e4));
%     %     i = i + 1;
%
%     %     R = GetBurst(R);
%     %     H = [H R.LFP.GammaBurstEvent.burst_du_steps{:}];
%     %     L = [L R.LFP.GammaBurstEvent.flat_du_steps{:}];
%
%     %     fprintf('\t File name: %s\n', files2{1}); % id_out-11
%     %     load(files2{1},'StiNeu')
%     switch ceil(id_out/100)
%         case 1
%             StiNeu = StiNeu1;
%         case 2
%             StiNeu = StiNeu2;
%         case 3
%             StiNeu = StiNeu3;
%         case 4
%             StiNeu = StiNeu4;
%     end
%     %%
%     switch 7 %R.ExplVar.NumP
%         case 1
%             Coor = [0;0];
%         case 2
%             Coor = [-16 15.5;-16 15.5];
%         case 3
%             Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
%         case 4
%             Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
%         case 5
%             Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
%         case 6
%             Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%                 -19.6 -9.8         -9.8        9.8          9.8         19.6];
%         case 7
%             Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%                 -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
%     end
%     hw = 31;
%     [Lattice, ~] = lattice_nD(2, hw);
%     LoalNeu = cell(1,R.ExplVar.NumP);
%     for i = 1:R.ExplVar.NumP
%         dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
%         LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
%     end
%     RasterPlotYL2(R,LoalNeu)
%     %%
%     RasterPlotYL2(R,StiNeu)
%
%     %     FR = vec2mat(R.num_spikes{1},10);
%     %     FR = sum(FR,2)'/3969*1e3;
%     %     FRcertain = filter(b,a,FR);
%     %     gamma_hilbert_abs = abs(hilbert(FRcertain)); % conv(abs(hilbert(FRcertain)), Kernel,'same');
%     %     GammaBurstEvent = GetBurstFR(gamma_hilbert_abs);
%     %     H = [H GammaBurstEvent.burst_du_steps{1}];
%     %     L = [L GammaBurstEvent.flat_du_steps{1}];
%
%     %     R = get_grid_firing_centreYL(R);
%     %     R = get_grid_firing_centreYL(R,'mode','bayesian');
%     %     R = {R};
%     %     SaveRYG(R);
%
%     %     IndC = [513 2497];
%     %     duration1 = 0.5e3:4.5e3;
%     %     duration2 = 5.5e3:9.5e3;
%     %     mode = 'quick';
%     %     subplot(2,3,1)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration1)
%     %     subplot(2,3,4)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration2)
%     %     mode = 'bayesian';
%     %     subplot(2,3,2)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration1)
%     %     subplot(2,3,5)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration2)
%     %     mode = 'bayesian2';
%     %     subplot(2,3,3)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration1)
%     %     subplot(2,3,6)
%     %     VisualizeGridFiringCentre2(R,mode,IndC,duration2)
%     %     next = input('\t Next figure?');
%     %     close all
%
%     saveas(gcf,[sprintf('%04g', id_out),'WMRasterPlot.eps'])
%     disp('Done');
%
%     %     bayes = R.grid.bayes;
%     %     save([sprintf('%04g-', id_out),'BWin10ms.mat'],'bayes')
% end
% Cenxy based on different items
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
repeat = 1000;
Cenx = cell(7,repeat);
Ceny = cell(7,repeat);
for NumP = 1:7
    switch NumP
        case 1
            Coor = [0;0];
        case 2
            Coor = [-16 15.5;-16 15.5];
        case 3
            Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
        case 4
            Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
        case 5
            Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
        case 6
            Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
                -19.6 -9.8         -9.8        9.8          9.8         19.6];
        case 7
            Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
                -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
    end
    LoalNeu = cell(1,NumP);
    R = load(files{1+repeat*(NumP-1)});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for id_out = 1:repeat
        fprintf('Processing output file No.%d out of %d...\n', id_out+repeat*(NumP-1), num_files);
        fprintf('\t File name: %s\n', files{id_out+repeat*(NumP-1)});
        R = load(files{id_out+repeat*(NumP-1)});
%         cenx = cell(1,NumP);
%         ceny = cell(1,NumP);
        for no = 1:NumP
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
            candx = Lattice(LoalNeu{no},1);
            candy = Lattice(LoalNeu{no},2);
            for i = 2.26e4:length(r)-bin/2 % ind % 4.46e4 % length(r) %
                if r(i) >= 120 % floor(0.5*length(candx)) %
                    % neurons within WM area
                    spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
%                     cenx{no} = [cenx{no} x-Coor(1,no)];
%                     ceny{no} = [ceny{no} y-Coor(2,no)];
                    
Cenx{NumP,id_out} = [Cenx{NumP,id_out} x-Coor(1,no)];
Ceny{NumP,id_out} = [Ceny{NumP,id_out} y-Coor(2,no)];
                end
            end
        end
%         Cenx{NumP,id_out} = cenx;
%         Ceny{NumP,id_out} = ceny;
    end
end
save('CenxyAll2.mat','Cenx','Ceny')
%% Cenxy based on different IE ratio??
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
repeat = 100;
Cenx = cell(10,repeat);
Ceny = cell(10,repeat);
NumP = 2;
for trial = 1:10
    Coor = [-16 15.5;-16 15.5];
    Coor(:,2) = Coor(:,2)-3*trial+2;
    LoalNeu = cell(1,NumP);
    R = load(files{1+100*(trial-1)});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for id_out = 1:repeat
        fprintf('Processing output file No.%d out of %d...\n', id_out+repeat*(trial-1), num_files);
        fprintf('\t File name: %s\n', files{id_out+repeat*(trial-1)});
        R = load(files{id_out+repeat*(trial-1)});
%         cenx = cell(1,NumP);
%         ceny = cell(1,NumP);
        for no = 1:NumP
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
            candx = Lattice(LoalNeu{no},1);
            candy = Lattice(LoalNeu{no},2);
            for i = 2.26e4:length(r)-bin/2 % ind % 4.46e4 % length(r) %
                if r(i) >= 120 % floor(0.5*length(candx)) %
                    % neurons within WM area
                    spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
%                     cenx{no} = [cenx{no} x-Coor(1,no)];
%                     ceny{no} = [ceny{no} y-Coor(2,no)];
                    
Cenx{trial,id_out} = [Cenx{trial,id_out} x-Coor(1,no)];
Ceny{trial,id_out} = [Ceny{trial,id_out} y-Coor(2,no)];
                end
            end
        end
%         Cenx{NumP,id_out} = cenx;
%         Ceny{NumP,id_out} = ceny;
    end
end
save('CenxyAll2.mat','Cenx','Ceny')
%% Cenxy based on different distances
bin = 500; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
repeat = 100;
Cenx = cell(10,repeat);
Ceny = cell(10,repeat);
NumP = 2;
for trial = 1:10
    Coor = [-16 15.5;-16 15.5];
    Coor(:,2) = Coor(:,2)-3*trial+2;
    LoalNeu = cell(1,NumP);
    R = load(files{1+100*(trial-1)});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for id_out = 1:repeat
        fprintf('Processing output file No.%d out of %d...\n', id_out+repeat*(trial-1), num_files);
        fprintf('\t File name: %s\n', files{id_out+repeat*(trial-1)});
        R = load(files{id_out+repeat*(trial-1)});
%         cenx = cell(1,NumP);
%         ceny = cell(1,NumP);
        for no = 1:NumP
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},:)),bin,2));
            candx = Lattice(LoalNeu{no},1);
            candy = Lattice(LoalNeu{no},2);
            for i = 2.26e4:length(r)-bin/2 % ind % 4.46e4 % length(r) %
                if r(i) >= 120 % floor(0.5*length(candx)) %
                    % neurons within WM area
                    spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2);
                    spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(LoalNeu{no},i-bin/2+1:i+bin/2));
                    [x,y,~,~,~,~] = fit_bayesian_bump_2_spikes_circular(spike_x_pos_o,spike_y_pos_o,2*hw+1,'quick');
%                     cenx{no} = [cenx{no} x-Coor(1,no)];
%                     ceny{no} = [ceny{no} y-Coor(2,no)];
                    
Cenx{trial,id_out} = [Cenx{trial,id_out} x-Coor(1,no)];
Ceny{trial,id_out} = [Ceny{trial,id_out} y-Coor(2,no)];
                end
            end
        end
%         Cenx{NumP,id_out} = cenx;
%         Ceny{NumP,id_out} = ceny;
    end
end
save('CenxyAll2.mat','Cenx','Ceny')
end

