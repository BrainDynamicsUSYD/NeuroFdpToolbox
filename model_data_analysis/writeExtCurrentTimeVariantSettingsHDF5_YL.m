function writeExtCurrentTimeVariantSettingsHDF5_YL(FID,pop_ind,mean,std,mean_TV,std_TV, meanTWO,mean_TVTWO,meanTHREE,mean_TVTHREE)
% adapt from function writeExtCurrentTimeVariantSettingsHDF5.m
% write external current settings with time-variant factors
%       FID: file id for writing data
%   pop_ind: neuron population index
%      mean: mean value for Gaussian currents (nA) for each neuron
%       std: std for Gaussrian currents for each neuron
%   mean_TV: a time-variant factor multiplied to the mean value for
%   Gaussian currents (nA) for each neuron (matrix)
%  mean_std: a time-variant factor multiplied to the std value for
%  Gaussrian currents for each neuron
%
% Explanation: For neuron i, at time step t, the external current will be a
% Gaussian random number with the sum(now:three) of a mean of mean[i]*mean_TV[t] and a std of std[i]*std_TV[t]

if length(mean) == 1 || length(std) == 1
    warning('INIT004 has been updated. MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;

h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean'],length(mean));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/std'],length(std));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TV'],length(mean_TV));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/std_TV'],length(std_TV));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/meanTWO'],length(meanTWO));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TVTWO'],length(mean_TVTWO));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/meanTHREE'],length(meanTHREE));
h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TVTHREE'],length(mean_TVTHREE));

h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean'],mean);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/std'],std);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TV'],mean_TV);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/std_TV'],std_TV);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/meanTWO'],meanTWO);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TVTWO'],mean_TVTWO);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/meanTHREE'],meanTHREE);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT019/mean_TVTHREE'],mean_TVTHREE);
end
