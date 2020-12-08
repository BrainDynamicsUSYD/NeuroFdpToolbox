function writeSPHDF5(FID,pre_pop_ind,post_pop_ind,SP_step,tau_D,tau_F,U)
% add SP to connection
%     FID: file id for writing data

pre_pop_ind = pre_pop_ind - 1; % from matlab to c++ index
post_pop_ind = post_pop_ind - 1;

syn_type = 0;  % AMPA
 
 
n_syns = h5read(FID,'/config/syns/n_syns');
n_match = NaN; 
for n = (1:n_syns)-1  % c++ index
    try
        type = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/type']);
        i_pre = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/i_pre']);
        j_post = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/j_post']);
        
        if pre_pop_ind == i_pre && post_pop_ind == j_post &&  syn_type == type
            n_match = n;
        end
    catch ME
    end
end


if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post and syn_type == AMPA!')
else
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT015/SP_on_step'], SP_step,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT015/tau_D'], tau_D,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT015/tau_F'], tau_F,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT015/U'], U,'WriteMode','append');
end

end

