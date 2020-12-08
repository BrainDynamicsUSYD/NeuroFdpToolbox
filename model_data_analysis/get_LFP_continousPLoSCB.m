function get_LFP_continousPLoSCB(varargin)
% adjust from function get_LFP_continousYL3.m
% Using LFP proxy on reference: Mazzoni, Alberto, et al. "Computing the
%      local field potential (LFP) from integrate-and-fire network models."
%      PLoS computational biology 11.12 (2015): e1004584.
dir_strut = dir('*0_neurosamp.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
fw = 63;
hw = 31;
LFP_range_sigma = 8;
% gamma power
fs = 1e4; % sampling frequency (Hz)
% Butterworth filter
order = 4; % 4th order
lowFreq = 30; % gamma band (default values for this function are 150-250 Hz)
hiFreq = 80;
Wn = [lowFreq  hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
% Loop number for PBS array job
loop_num = 0;
%%
% for k = [2 11:num_files]% 1:num_files
%     loop_num = loop_num + 1;
%     % For PBS array job
%     if nargin ~= 0
%         PBS_ARRAYID = varargin{1};
%         if loop_num ~=  PBS_ARRAYID
%             continue;
%         end
%     end
%     R = load(files{k});
%     disp('Creating LFP_grid...')
%     [~, steps] = size(R.I_AMPA);
%     [Lattice, ~] = lattice_nD(2, hw);
%     dist = lattice_nD_find_dist(Lattice, hw,  0, 0);
%     gaus_tmp = 1/(LFP_range_sigma*sqrt(2*pi))*exp(-0.5*(dist/LFP_range_sigma).^2) .* double(dist <= LFP_range_sigma*2.5);
%     m = reshape(gaus_tmp, fw, fw);
%     LFP_grid = zeros(fw,fw,steps);
%     for i = 61:steps
%         x = reshape(abs(R.I_AMPA(:,i-60)) + 1.65*abs(R.I_GABA(:,i)), fw, fw);
%         y = convolve2(x, m, 'wrap');
%         LFP_grid(:,:,i) = y;
%     end
%     LFP_grid = LFP_grid(:,:,61:end);
%     gamma_power_gridPLoS = zeros(size(LFP_grid));
%     disp('Creating gamma_power_gridPLoS...')
%     for i = 1:fw
%         for j= 1:fw
%             gamma_tmp = filter(b,a,LFP_grid(i,j,:));
%             gamma_tmp = gamma_tmp(:);
%             hil_tmp = abs(hilbert( gamma_tmp));
%             gamma_power_gridPLoS(i,j,:) = hil_tmp;
%         end
%     end
%     clear LFP_grid;
%     save(files{k}, 'gamma_power_gridPLoS', '-append')
% end
%%
[Lattice, ~] = lattice_nD(2, hw);
S1 = 3.1; % 20 for default
Local1 = {};
LFP_centre_x = linspace(-hw, hw, 21); % E16:9  E100:21 E400:41
LFP_centre_y = linspace(-hw, hw, 21);
LFP_centre_x = LFP_centre_x(2:2:20); % E16(2:2:8)  E100(2:2:20) E400(2:2:40)
LFP_centre_y = LFP_centre_y(2:2:20);
[LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
LFP_centre_x = LFP_centre_x(:);
LFP_centre_y = LFP_centre_y(:);
for i = 1:length(LFP_centre_x)
    %     Local1{i} = find(lattice_nD_find_dist(Lattice,hw,LFP_centre_x(i),LFP_centre_y(i)) <= S1);
    x1 = LFP_centre_x(i) - S1;
    x2 = LFP_centre_x(i) + S1;
    y1 = LFP_centre_y(i) - S1;
    y2 = LFP_centre_y(i) + S1;
    Local1{i} = find(Lattice(:,1)>=x1 & Lattice(:,1)<=x2 & Lattice(:,2)>=y1 & Lattice(:,2)<=y2);
end
for k = [2 11:num_files]% 1:num_files
    loop_num = loop_num + 1;
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    load(files{k},'I_AMPA','I_GABA');
    disp('Creating LFP_grid...')
    [~, steps] = size(I_AMPA);    
    LFP_temp = zeros(100,steps);
    for i = 61:steps
        for j = 1:length(LFP_centre_x)
            LFP_temp(j,i) = sum(abs(I_AMPA(Local1{j}(:),i-60)) + 1.65*abs(I_GABA(Local1{j}(:),i)));
        end
    end
    LFP_grid = reshape(LFP_temp,[10 10 steps]);
    LFP_grid = LFP_grid(:,:,61:end);
    gamma_power_grid = zeros(size(LFP_grid));
    disp('Creating gamma_power_gridPLoS...')
    for i = 1:10
        for j= 1:10
            gamma_tmp = filter(b,a,LFP_grid(i,j,:));
            gamma_tmp = gamma_tmp(:);
            hil_tmp = abs(hilbert( gamma_tmp));
            gamma_power_grid(i,j,:) = hil_tmp;
        end
    end
    clear LFP_grid;
    save([sprintf('%04g-',k),'ModifitedSquareLFPGammaGrid.mat'],'LFP_temp','gamma_power_grid')
%     save(files{k}, 'gamma_power_grid2', '-append')
end
end