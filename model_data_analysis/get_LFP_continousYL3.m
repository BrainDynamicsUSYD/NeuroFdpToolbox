function get_LFP_continousYL3(varargin)
% adjust from function get_LFP_continousYL2.m
% only work for gamma_power_grid/gamma_phase_grid
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
lowFreq = 100; % gamma band (default values for this function are 150-250 Hz)
hiFreq = 250;
Wn = [lowFreq  hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
% Loop number for PBS array job
loop_num = 0;
for k = 1:num_files % 1:num_files
    loop_num = loop_num + 1;
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    R = load(files{k},'I_AMPA','I_GABA');
    disp('Creating LFP_grid...')
    [~, steps] = size(R.I_AMPA);
    [Lattice, ~] = lattice_nD(2, hw);
    dist = lattice_nD_find_dist(Lattice, hw,  0, 0);
    gaus_tmp = 1/(LFP_range_sigma*sqrt(2*pi))*exp(-0.5*(dist/LFP_range_sigma).^2) .* double(dist <= LFP_range_sigma*2.5);
    m = reshape(gaus_tmp, fw, fw);
    LFP_grid = zeros(fw,fw,steps);
    for i = 1:steps
        x = reshape(abs(R.I_AMPA(:,i)) + abs(R.I_GABA(:,i)), fw, fw);
        y = convolve2(x, m, 'wrap');
        LFP_grid(:,:,i) = y;
    end
    clear x y R
    gamma_power_grid = zeros(size(LFP_grid));
    disp('Creating gamma_power_grid...')
    for i = 1:fw
        for j= 1:fw
            gamma_tmp = filter(b,a,LFP_grid(i,j,:));
            gamma_tmp = gamma_tmp(:);
            gamma_power_grid(i,j,:) = abs(hilbert( gamma_tmp));
%             gamma_phase_grid(i,j,:) = angle(hilbert( gamma_tmp));
        end
    end
    SWR_power_grid = gamma_power_grid;
    clear LFP_grid gamma_power_grid
    save(files{k}, 'SWR_power_grid','-append') % ,'-v7.3') %, '-append')
%     save('0006-GammaPhaseGrid.mat', 'gamma_phase_grid','-v7.3')
end
end