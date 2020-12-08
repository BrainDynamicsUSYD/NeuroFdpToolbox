function main_toymodel1()
% Use coherent units (msec+mV+nF+miuS+nA) unless otherwise stated
% using matlab only to test the equations
%%%% Seed the Matlab random number generator
seed = 1;
rng(seed);

% generate time string with a precision up to msec
date = datestr(now,'dd-mm-yyyy HH:MM:SS FFF');

%%%% Define some basic parameters
% Time step (ms)
dt = 0.1; % ms
% Total number of simulation steps
step_tot = 10^5; % 10 s
N = 2;

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
V_old1 = V_rt;
V_old2 = V_rt;
Ref_steps_left1 = 0;
Ref_steps_left2 = 0;
g_K1 = 0;
g_K2 = 0;
tau_K = 80; % ms
dg_K = 0.01; % miuS
Vec_spikes1 = false(1,step_tot);
Vec_spikes2 = false(1,step_tot);

%%%% Define the initial condition
p_fire = [0.1]; %  0.1]; % initial firing probabilities for both populations
% set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
r_V0 = [0.5]; %  0.5];

%%%% Add external Poissonian spike trains into the 1st population
pop = 1;
type_ext = 1; % 1 for AMPA-like syanpses
g_ext = 2*10^-3; % synaptic coupling strength (miuS)
N_ext = 1000; % No. of independent external connections onto each neuron
rate_ext = 0.6;
neurons_recv = zeros(1, N(pop));
neurons_recv(1) = 1; % only neurons #1-10 receive such external input

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

%%%% Add synaptic plasticity to the synapses with the 1st population
pop_pre = 1;
pop_post = 1;
step_start = 5e3; % round(step_tot/3); % Turn on SP at this time step


% simulate
wb_h = waitbar(0,'Please wait....');
for t = 1:step_tot-1
    g_K1 = (1-dt/tau_K)*g_K1;
    g_K2 = (1-dt/tau_K)*g_K2;
    if Ref_steps_left1 > 0
        Ref_steps_left1 = Ref_steps_left1 - 1;
    end
    if Ref_steps_left2 > 0
        Ref_steps_left2 = Ref_steps_left2 - 1;
    end    
    g_K1 = g_K1 + dg_K*double(Vec_spikes1(t));
    g_K2 = g_K2 + dg_K*double(Vec_spikes2(t));
    if mod(t,round(step_tot/100)) == 0
       waitbar(t/step_tot, wb_h);
    end
    AJM = AJ*double(Mat_spikes(:,t));% matlab-specific operation
    ext_spike = poissrnd(rate_ext*N_ext*dt, 1, 1);
    V_new1 =  (1-g_lk*dt/Cm)*V_old1 + g_lk*dt/Cm*V_lk + ext_spike*j  + V_old;% matlab-specific operation: matrix operations instead nested for-loops
    V_new(Ref_steps_left1 >  0) = V_r; % matlab-specific operation: logical indexing
    Mat_spikes(V_new >= V_thr,t+D/dt) = true;
    Ref_steps_left1(V_new >= V_thr) = round(t_rp/dt);
%     V_sample(:,t+1) = V_new(i_sample);
    V_old = V_new;
%     [I_exc,J_exc] = ind2sub(s,find(Mat_spikes(1:N_E,t)));
%     h = plot(I_exc,J_exc,'bo');
%     axis([0 l 0 l]);
%     pause(0.01);
%     delete(h);
end
close(wb_h);

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