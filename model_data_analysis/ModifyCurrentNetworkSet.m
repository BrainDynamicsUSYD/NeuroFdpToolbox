function ModifyCurrentNetworkSet(varargin)
PBS_ARRAYID = varargin{1};
dir_strut = dir([sprintf('%04g-', PBS_ARRAYID),'*in.h5']);
FID = dir_strut.name;
dir_strut2 = dir('*_config_data.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end
load(files2{1},'StiNeu');

hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=44 gives 7921
N_e = (hw*2+1)^2;
Mean = 3*ones(1,N_e); % 5
Std = zeros(1,N_e);
NumP = 4;
dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 10*sec; % use 120 seconds!
mean_TV = zeros(NumP+1,step_tot);
std_TV = zeros(NumP+1,step_tot);
TV_group = (NumP+1)*ones(1,N_e);
start = 2e4; % 0.1 ms
mean_TV(PBS_ARRAYID,start:(start+0.25e4)) = 1;
TV_group(StiNeu{PBS_ARRAYID}) = PBS_ARRAYID;
writeExtCurrentTimeVariantSettingsMultiGroupHDF5_subfunction(FID, 1, Mean, Std, mean_TV, std_TV, TV_group)
disp('Finishing overwriting for loading item.')
end