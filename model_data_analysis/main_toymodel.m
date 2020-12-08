function main_toymodel()
% Use coherent units (msec+mV+nF+miuS+nA) unless otherwise stated

%%%% Seed the Matlab random number generator
% if nargin ~= 0
%     PBS_ARRAYID = varargin{1};
% end
seed = 3;

%%%% Open new (uniquely time-stamped) input files for the SpikeNet C++
% simulator.
[FID] = new_ygin_files_and_randseedHDF5(seed);

%%%% Define some basic parameters
% Time step (ms)
dt = 0.1; % 0.1 is usually small enough
% Total number of simulation steps
step_tot = 10^5;
% Create two neuron populations with 50 neurons in the 1st population
% and 10 neurons in the 2nd population.
N = [2];
% Write the above basic parameters to the input file
writeBasicParaHDF5(FID, dt, step_tot, N);

%%%% Define non-parameters for the neuorn model
% For demo purpose, the values used as following are the still the
% default ones.
% If default values are to be used, there is no need to re-define them.
Cm = 0.25; % (nF) membrane capacitance
tau_ref = 4.0; % (ms) absolute refractory time
V_rt = -60.0;   % (mV) reset potential
V_lk = -70.0;   % (mV) leaky reversal potential
V_th = -50.0;   % (mV) firing threshold
g_lk = 0.0167;   % (miuS) leaky conductance (note that Cm/gL=15 ms)
for pop_ind = 1 % :2
    writePopParaHDF5(FID, pop_ind,'Cm',Cm,'tau_ref',tau_ref,'V_rt',V_rt,...
        'V_lk',V_lk,'V_th',V_th,'g_lk',g_lk);
end

% %%%% Add spike-frequency adaptation to the 1st population
pop = 1;
writeSpikeFreqAdptHDF5(FID, pop);

% % Use Adam 2016 synapse model
% model_choice = 2;
% writeSynapseModelChoiceHDF5(FID, model_choice)

%%%% Define the initial condition
p_fire = [0.1]; %  0.1]; % initial firing probabilities for both populations
% set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0]
r_V0 = [0.5]; %  0.5];
writeInitCondHDF5(FID, r_V0, p_fire)

% %%%% Add external Poissonian spike trains into the 1st population
pop = 1;
type_ext = 1; % 1 for AMPA-like syanpses
g_ext = 2*10^-3; % synaptic coupling strength (miuS)
N_ext = 1000; % No. of independent external connections onto each neuron
rate_ext = 1*ones(1, step_tot); % Poisson rate for each time step (Hz)
rate_ext(5e4:7e4) = 1.5;
rate_ext(7e4:8e4) = 0.7;
neurons_recv = zeros(1, N(pop));
neurons_recv(1) = 1; % only neurons #1-10 receive such external input
writeExtSpikeSettingsHDF5(FID, pop, type_ext, g_ext,  N_ext, rate_ext,...
    neurons_recv );

% %%%% Add external currents (Gaussian white noise) to the 2nd population
% pop = 1;
% I_ext_mean = 0.5*[1 0]; % defined for each neuron (nA)
% I_ext_std = 0.5*[1 0]; % defined for each neuron (nA)
% mean_TV = ones(1,step_tot);
% mean_TV(5e4:7e4) = 1.5;
% mean_TV(7e4:8e4) = 0.7;
% std_TV = ones(1,step_tot);
% writeExtCurrentSettingsHDF5(FID, pop, I_ext_mean, I_ext_std)
% writeExtCurrentTimeVariantSettingsHDF5(FID, pop, I_ext_mean, I_ext_std, mean_TV, std_TV)
%
% %%%% Add external conductances (Gaussian white noise) to the 2nd population
% % The default reveral potential is V_ext = 0.0 mV.
% pop = 1;
% g_ext_mean = 50*1e-3*[1 0]; % defined for each neuron (muS)
% g_ext_std = 0.1*1e-3*[1 0]; % defined for each neuron (muS)
% writeExtConductanceSettingsHDF5(FID, pop, g_ext_mean, g_ext_std)
% writeExtConductanceTimeVariantSettingsHDF5(FID, pop, g_ext_mean, g_ext_std, mean_TV, std_TV)

%%%% Define runaway killer
% The computational cost of the simulation is directly proportional to
% the avearage firing rate of the network. Very often, your parameters
% may lead to biologically unrealistically high firing rate or even
% runaway activity (the highest rate possible ~1/tau_ref). To save the
% computational resources for more interesting simulation cases, it is
% advised to define the runawary killer.
% Note that the simulator will still output all the relevant data before
% the kill.
min_ms = 500; % the min simulation duration that should be guaranteed
runaway_Hz = 40; % the threshold above which the simu should be killed
Hz_ms = 200; % the window length over which the firing rate is averaged
pop = 1; % the population to be monitored by the runaway killer
writeRunawayKillerHDF5(FID, pop, min_ms, runaway_Hz, Hz_ms);

%%%% Record the basic statistics for the 1st neuron population
pop = 1;
writePopStatsRecordHDF5(FID, pop);

%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMPA, 2:GABAa, 3:NMDA)
% Define AMPA-like excitatory synaptic connection wtihin the 1st pop
pop_pre = 1;
pop_post = 1;
syn_type = 1; % 1 for AMPA-like synapse
g_EE = 10*10^-3; % synaptic coupling strength (miuS)
% random adjacency matrix with a wiring probability of 0.2
% A = rand(N(1), N(1)) < 1; % 0.2;
% [I_EE, J_EE] = find(A);
I_EE = 1;
J_EE = 2;
K_EE = g_EE*ones(size(I_EE)); % identical coupling strength
D_EE = rand(size(I_EE))*4; % uniformly random conduction delay [0, 4] ms
writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
    I_EE, J_EE, K_EE, D_EE);

% % Define AMPA-like excitatory synaptic connection from pop 1 to pop 2
% pop_pre = 1;
% pop_post = 2;
% syn_type = 1; % 1 for AMPA-like synapse
g_IE = 10*10^-3; % synaptic coupling strength (miuS)
% % random adjaceny matrix with a wiring probability of 0.5
% A = rand(N(1), N(2)) < 0.5;
% [I_IE, J_IE] = find(A);
% K_IE = g_IE*ones(size(I_IE)); % identical coupling strength
% D_IE = rand(size(I_IE))*4; % uniformly random conduction delay [0, 4] ms
% writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
%     I_IE, J_IE, K_IE, D_IE);
%
% % Define GABA-like excitatory synaptic connection from pop 2 to pop 1
% pop_pre = 2;
% pop_post = 1;
% syn_type = 2; % 2 for GABA-like synapse
g_EI = 10*10^-3; % synaptic coupling strength (miuS)
% % random adjaceny matrix with a wiring probability of 0.5
% A = rand(N(2), N(1)) < 0.5;
% [I_EI, J_EI] = find(A);
% K_EI = g_EI*ones(size(I_EI)); % identical coupling strength
% D_EI = rand(size(I_EI))*4; % uniformly random conduction delay [0, 4] ms
% writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
%     I_EI, J_EI, K_EI, D_EI);
%
% % Define GABA-like excitatory synaptic connection within pop 2
% pop_pre = 2;
% pop_post = 2;
% syn_type = 2; % 2 for GABA-like synapse
g_II = 20*10^-3; % synaptic coupling strength (miuS)
% % random adjaceny matrix with a wiring probability of 0.5
% A = rand(N(2), N(2)) < 0.5;
% [I_II, J_II] = find(A);
% K_II = g_II*ones(size(I_II)); % identical coupling strength
% D_II = rand(size(I_II))*4; % uniformly random conduction delay [0, 4] ms
% writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
%     I_II, J_II, K_II, D_II);
%%%%%%%%%%%%%%%%%%% Chemical Connections Done %%%%%%%%%%%%%%%%%%%%%%%

%%%% Use an existing synapse connectivity definition file instead of
% generating a new one:
% (The feature needs to be written)

% %%%% Define non-default synapse parameters
% writeSynParaHDF5(FID, 'tau_decay_GABA', 6);
% % "help writeSynPara" to see all the parameter names and default values

% %%%% Add short-term depression to the synapses with the 1st population
% pop_pre = 1;
% pop_post = 1;
% step_start = round(step_tot/3); % Turn on STD at this time step
% writeSTDHDF5(FID, pop_pre, pop_post, step_start);

%%%% Add synaptic plasticity to the synapses with the 1st population
pop_pre = 1;
pop_post = 1;
step_start = 5e3; % round(step_tot/3); % Turn on STD at this time step
% writeSPHDF5(FID, pop_pre, pop_post, step_start);

% %%%% Add inhibitory STDP to the synapses from pop 2 to pop 1
% pop_pre = 2;
% pop_post = 1;
% % Turn on inhibitory STDP at a certain time step
% step_start = round(step_tot/3*2);
% writeInhSTDPHDF5(FID, pop_pre, pop_post, step_start);

%%%% Record the basic statistics for the excitatory synapses within
% the 1st neuron population
pop_pre = 1;
pop_post = 1;
syn_type = 1; % 1 for AMPA-like synapse
writeSynStatsRecordHDF5(FID, pop_pre, pop_post, syn_type);

%%%% Sample detailed time series data from the 1st population
pop = 1;
sample_neuron = 1:2; % 10:N(1); % sample every 10th neuron
sample_steps = 1:step_tot; % 2:step_tot; % sample every 2nd time step
sample_data_type = [1,1,1,1,0,0,1,1];
% The logical vector above corresponds to
% [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
writeNeuronSamplingHDF5(FID, pop, sample_data_type, ...
    sample_neuron, sample_steps)

%%%% Optional: record explanatory variables (scalars only)
% They will also be used in pro-processsing for auto-generated comments
discard_transient = 0; % transient period data to be discarded (ms)
comment1 = 'This is a demo.';
comment2 = ['If you have any question regarding SpikeNet,'...
    'please contact Yifan Gu (yigu8115@gmail.com).'];
writeExplVarHDF5(FID, 'discard_transient', discard_transient, ...
    'g_EE', g_EE, ...
    'g_EI', g_EI, ...
    'g_IE', g_IE, ...
    'g_II', g_II, ...
    'g_ext', g_ext, ...
    'N_ext', N_ext, ...
    'rate_ext', rate_ext,...
    'comment1', comment1, 'comment2', comment2);

%%%% append this matlab file to the input file for future reference
appendThisMatlabFileHDF5(FID)

end


function appendThisMatlabFileHDF5(FID)

% need to coopy and past the following code into the file to be appended!
text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,'/config/MATLAB/config.m',text,'WriteMode','append');

end