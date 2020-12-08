function [R] = get_CC_pop(R, pops)

fprintf('\t Getting correlation coefficient distribution for specified population(s)...\n');

corrcoef_sample_num = 10^4; % number of sampling pairs
CC_kernel_width = 40; % ms, kernel length

% Dumpe fields
N = R.N;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;

% Define kernel
kernel_type = 'square_Hz';
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);

% corrcoef within each population
CC_pop = cell(Num_pop,1);
CC_pop_pairs = cell(Num_pop,1);

if nargin == 1
    pops = 1:Num_pop;
end

for pop_ind = pops
    if nnz(num_spikes{pop_ind}) > 0

        pairs = rand_unique_pairs( N(pop_ind), corrcoef_sample_num );
        CC_pop_pairs{pop_ind} = pairs;
        % Calculate pair-wise corrcoef
        for i = 1:corrcoef_sample_num
            CC_pop{pop_ind}(1,i) = CorrCoefYG(spike_hist{pop_ind}(pairs(1,i),:), spike_hist{pop_ind}(pairs(2,i),:), CC_kernel);
            % display progress
            
%             if i > 1
%                 fprintf(repmat('\b',1,10));
%             end
%             if i < corrcoef_sample_num
%                 fprintf('%10g',i);
%             end
        end
     
    end
end

% Record results
R.Analysis.CC_pop_kernel_type = kernel_type;
R.Analysis.CC_pop_kernel_width = CC_kernel_width;
R.Analysis.CC_pop = CC_pop;
R.Analysis.CC_pop_pairs = CC_pop_pairs;

end


