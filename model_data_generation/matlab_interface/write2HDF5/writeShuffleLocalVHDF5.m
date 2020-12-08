function writeShuffleLocalVHDF5(FID, pop_ind, neuron_ind, time_index)
% write shuffle V of local area settings for reducing #patterns
%        FID: file id for writing data
%    pop_ind: neuron population index
% neuron_ind: matrix of neuron indices to be sampled in pop_ind, rows are
%             local area, columns are neurons.
% time_index: logical vector to define time points to be sampled
%
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1,0]
% 4,000 neurons x 10,000 time points will be 300MB data


% if sum((time_index ~= 0) & (time_index ~= 1)) ~= 0
%     error('sample_time_index must be logical vectors: ');
% end

% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SHUFFLE_V/neurons'],neuron_ind,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SHUFFLE_V/time_points'],time_index,'WriteMode','append');
end
