function main_Gamma_SP2(varargin)
% adjust from function main_SWR_reference_gaus_balance.m

dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 10*sec;
discard_transient = 0; % ms

% Loop number for PBS array job
loop_num = 0;

tau_ref = 4; % ms
delay = 4; % ms

% sptially embedded network
hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=44 gives 7921
N_e = (hw*2+1)^2;
N_i = 1000;
N = [N_e, N_i];
Num_pop = length(N);
Type_mat = ones(Num_pop);
Type_mat(end, :) = 2;
in_out_r = [0.13];
% [Lattice,~] = lattice_nD(2, hw);

% synaptic plasticity
STD_on = 0;
SP_on = 1; % need to get two subfunctions consistent

repeats = 100;

for Kw = [1.8] % 2
    for ring = [20]
        for NumP = [1:4]           
            for ring_wid = [8.2]
                switch NumP
                    case 1
                        Coor = [0;-ring];
                    case 2
                        Coor = [0 0;-ring ring];
                    case 3
                        Coor = [-ring*sqrt(3)/2 ring*sqrt(3)/2 0;-ring/2 -ring/2 ring];
                    case 4
                        Coor = [0 -ring ring 0;-ring 0 0 ring];
                    case 8
                        Coor = [0 -ring/sqrt(2) ring/sqrt(2) -ring ring -ring/sqrt(2) ring/sqrt(2) 0;...
                            -ring -ring/sqrt(2) -ring/sqrt(2)  0     0   ring/sqrt(2) ring/sqrt(2) ring];
                end
                Stogether = 1;
                % randsample(N_e,NumP); overlap[513 2497 2503]
                mean_TV = zeros(NumP+1,step_tot);
                std_TV = zeros(NumP+1,step_tot);
                Std = zeros(1,N_e);
                TV_group = (NumP+1)*ones(1,N_e);
                tau_D = 200; % 200 ms
                tau_F = 1500; % 1500 ms
                U = 0.2;
                for P0_init = 0.08*ones(1,repeats)
                    P_mat = [P0_init 0.1;
                        0.1  0.2]*2;
                    
                    % parameter
                    for decay_GABA = [6] % 6
                        for decay_AMPA = [5] % 5
                            for SpikeFreqAapt = [0]
                                for tau_K = [80] % ms 80
                                    for dg_K = 0.01 % uS 0.01
                                        for LFP_range_sigma = [8] % 8
                                            for cn_scale_wire = [2] % 2
                                                for cn_scale_weight = [2] % 2
                                                    iter_num = 5;
                                                    N_ext = 1000;
                                                    g_ext = 2*10^-3; % uS
                                                    [ fit_g_2_EPSP_2, ~ ] = g_EPSP_conversion( );
                                                    for deg_hybrid = [0.4]
                                                        degree_CV = 0.2; % 0.2 works
                                                        for g_mu = [5]*10^-3 % 5 4 15
                                                            EPSP_mu = fit_g_2_EPSP_2(g_mu);
                                                            EPSP_sigma = 1;
                                                            for inh_STDP = [0]
                                                                
                                                                %  K_ee_mean is about 0.5, need 1000 in-coming connections.
                                                                %  this is not good.
                                                                %  what can I do???
                                                                %  ref: A Lognormal Recurrent Network Model for Burst Generation during Hippocampal Sharp Waves
                                                                
                                                                for g_balance = [1] % 1
                                                                    for g_EI = [ 13.5 ]*10^-3 % 13.5
                                                                        for g_IE = [5]*10^-3 % 5
                                                                            for g_II = [25]*10^-3 % 25
                                                                                for rate_ext_I = [1] % 1 kHz
                                                                                    for rate_ext_E = [0.85] % 0.85
                                                                                        for tau_c_EE = [8] % the three are exponential distance dependent time constants
                                                                                            tau_c_IE = 10;
                                                                                            for tau_c_I = [20]
                                                                                                for backgroundcoe = [1]
                                                                                                    backgroundE = [1]; % 1:normal,2:modify
                                                                                                    
                                                                                                    % For PBS array job
                                                                                                    loop_num = loop_num + 1;
                                                                                                    if nargin ~= 0
                                                                                                        PBS_ARRAYID = varargin{1};
                                                                                                        if loop_num ~=  PBS_ARRAYID
                                                                                                            continue;
                                                                                                        end
                                                                                                    end
                                                                                                    
                                                                                                    % seed the matlab rand function! The seed is global.
                                                                                                    [FID] = new_ygin_files_and_randseedHDF5(loop_num);
                                                                                                    
                                                                                                    % write basic parameters
                                                                                                    writeBasicParaHDF5(FID, dt, step_tot, N);
                                                                                                    
                                                                                                    K_mat = [NaN  g_IE;
                                                                                                        g_EI  g_II]; % miuSiemens
                                                                                                    
                                                                                                    if SpikeFreqAapt == 1
                                                                                                        writeSpikeFreqAdptHDF5(FID, 1,dg_K,tau_K);
                                                                                                    end
                                                                                                    
                                                                                                    % write pop para
                                                                                                    for pop_ind = 1:Num_pop
                                                                                                        writePopParaHDF5(FID,pop_ind,'tau_ref',tau_ref);
                                                                                                    end
                                                                                                    
                                                                                                    % write external currents
                                                                                                    % the third paramster value 1 means that external inputs are AMPA currensts
                                                                                                    background = rate_ext_E*ones(1, step_tot);
                                                                                                    if backgroundE == 2
                                                                                                        background(2e4+1:end) = backgroundcoe*background(2e4+1:end);
                                                                                                    end
                                                                                                    writeExtSpikeSettingsHDF5(FID, 1, 1, g_ext,  N_ext, background,  ones(1, N(1)) );
                                                                                                    writeExtSpikeSettingsHDF5(FID, 2, 1, g_ext,  N_ext, rate_ext_I*ones(1, step_tot),  ones(1, N(2)) );
                                                                                                    
                                                                                                    % write synapse para
                                                                                                    writeSynParaHDF5(FID,'tau_decay_GABA',decay_GABA,'tau_decay_AMPA',decay_AMPA); % SWR:3;5
                                                                                                    
                                                                                                    % inhibitory STDP
                                                                                                    if inh_STDP == 1
                                                                                                        writeInhSTDPHDF5(FID, 2, 1, 1*sec);
                                                                                                    end
                                                                                                    
                                                                                                    if STD_on == 1
                                                                                                        writeSTDHDF5(FID, 1, 1, 1);
                                                                                                    end
                                                                                                    
                                                                                                    % write runaway killer
                                                                                                    min_ms = 500; % 5 sec
                                                                                                    runaway_Hz = 100; % ??
                                                                                                    Hz_ms = 200; % ms
                                                                                                    writeRunawayKillerHDF5(FID, 1, min_ms, runaway_Hz, Hz_ms);
                                                                                                    
                                                                                                    % random initial condition settings (int pop_ind, double p_fire)
                                                                                                    p_fire = [0.1 0.00]; % between [0,1], 0.05
                                                                                                    writeInitVHDF5(FID, p_fire);
                                                                                                    
                                                                                                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                                                                                    % type(1:AMPA, 2:GABAa, 3:NMDA)
                                                                                                    
                                                                                                    Lattice_I = quasi_lattice_2D( N(2) , hw);
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                                                                    % generate in- and out-degree accroding to (1) distance-dependent rule
                                                                                                    % (2) degree distributios and (3) common neighbour rule
                                                                                                    deg_mean = N_e*P0_init;
                                                                                                    deg_std_logn = degree_CV*deg_mean;
                                                                                                    
                                                                                                    [deg_in_0,deg_out_0] = hybrid_degree(N_e,deg_mean,deg_std_logn,in_out_r,deg_hybrid);
                                                                                                    
                                                                                                    [I_ee,J_ee,dist_IJ,iter_hist,Lattice_E] = generate_IJ_2D(deg_in_0,deg_out_0,tau_c_EE,cn_scale_wire,iter_num);
                                                                                                    in_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 1)); % row vector
                                                                                                    out_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 2));
                                                                                                    
                                                                                                    % Generate K according to
                                                                                                    mu_p = log((EPSP_mu^2)/sqrt(EPSP_sigma^2+EPSP_mu^2));
                                                                                                    s_p = sqrt(log(EPSP_sigma^2/(EPSP_mu^2)+1));
                                                                                                    mu_p = mu_p + s_p^2;
                                                                                                    g_pool_generator_hld = @(N)g_pool_generator(N,mu_p,s_p);
                                                                                                    K_scale = sqrt(in_degree);
                                                                                                    K_cell = inverse_pool( in_degree, K_scale, g_pool_generator_hld);
                                                                                                    K_ee = NaN;
                                                                                                    if ~isnan(K_cell{1})
                                                                                                        K_ee = zeros(size(J_ee));
                                                                                                        for j = 1:N_e
                                                                                                            K_ee(J_ee==j) = K_cell{j}';
                                                                                                        end
                                                                                                        clear K_cell; % reformat K
                                                                                                    end
                                                                                                    
                                                                                                    % shuffle K accordind to common neighbour rule
                                                                                                    if ~isnan( K_ee)
                                                                                                        [K_ee] = shuffle_K_common_neighbour(K_ee,I_ee,J_ee,cn_scale_weight);
                                                                                                    end
                                                                                                    %%% change EE coupling for local populations
                                                                                                    [K_ee,RingNeu] = WorkingMemoryLocalPopulation2(hw,Lattice_E,I_ee,J_ee,K_ee,Kw,ring,ring_wid);
                                                                                                    K_ee_mean = mean(K_ee);
                                                                                                    EE_input = full(sum(sparse(I_ee,J_ee,K_ee),1));
                                                                                                    
                                                                                                    D_ee = rand(size(I_ee))*delay;
                                                                                                    writeChemicalConnectionHDF5(FID,Type_mat(1,1),1,1,I_ee,J_ee,K_ee,D_ee);
                                                                                                    clear I J K D;
                                                                                                    
                                                                                                    % [~,ind_sorted] = sort(in_degree);
                                                                                                    % sample_neuron = ind_sorted(1:500:end);
                                                                                                    % sample_neurons = [StiNeu{1}(1:5)' StiNeu{2}(1:5)' StiNeu{3}(1:5)' 1:5];
                                                                                                    StiNeu = cell(1,NumP);
                                                                                                    for i = 1:NumP
                                                                                                        dist = Distance_xy(Lattice_E(RingNeu,1),Lattice_E(RingNeu,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
                                                                                                        StiNeu{i} = RingNeu(dist<=ring_wid/2)';
                                                                                                    end
                                                                                                    OutNeu = setdiff(1:N_e,[StiNeu{:}]); % row vector
                                                                                                    sample_neurons = [];
                                                                                                    for i = 1:length(StiNeu)
                                                                                                        sample_neurons = [sample_neurons StiNeu{i}(1)']; % StiNeu{i} column vector
                                                                                                    end
                                                                                                    size(sample_neurons)
                                                                                                    sample_neurons = [sample_neurons OutNeu(1)];
                                                                                                    
                                                                                                    Mean = 3*ones(1,N_e); % 5
                                                                                                    start = 2e4+1; % 0.1 ms
                                                                                                    for i = 1:NumP
                                                                                                        mean_TV(i,start:(start+0.25e4)) = 1;
                                                                                                        if Stogether == 0
                                                                                                            start = start + 2e4; % 2e4
                                                                                                        end
                                                                                                        TV_group(StiNeu{i}) = i;
                                                                                                    end
                                                                                                    writeExtCurrentTimeVariantSettingsMultiGroupHDF5(FID, 1, Mean, Std, mean_TV, std_TV, TV_group)
                                                                                                    
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                                                                    [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_E, hw, tau_c_I, P_mat(2,1) );
                                                                                                    D = rand(size(I))*delay;
                                                                                                    K = zeros(size(J));
                                                                                                    for i_E = 1:N(1)
                                                                                                        mu_K_tmp = EE_input(i_E)/sum(J==i_E)*(g_EI/g_mu)*g_balance; % g_mu
                                                                                                        K(J==i_E) = abs(randn([1 sum(J==i_E)])*(mu_K_tmp/4) + mu_K_tmp); % this is a bit too arbitary!
                                                                                                    end
                                                                                                    writeChemicalConnectionHDF5(FID, Type_mat(2, 1),  2, 1,   I,J,K,D); % EI
                                                                                                    clear I J K D;
                                                                                                    
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                                                                    [ I,J ] = Lattice2Lattice( Lattice_E, Lattice_I, hw, tau_c_IE, P_mat(1,2) );
                                                                                                    D = rand(size(I))*delay;
                                                                                                    K = ones(size(I))*K_mat(1,2);
                                                                                                    writeChemicalConnectionHDF5(FID, Type_mat(1, 2),  1, 2,   I,J,K,D);
                                                                                                    clear I J K D;
                                                                                                    
                                                                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                                                                    [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_I, hw, tau_c_I, P_mat(2,2) );
                                                                                                    D = rand(size(I))*delay;
                                                                                                    K = ones(size(I))*K_mat(2,2);
                                                                                                    writeChemicalConnectionHDF5(FID, Type_mat(2, 2),  2, 2,   I,J,K,D);
                                                                                                    clear I J K D;
                                                                                                    
                                                                                                    if SP_on == 1
                                                                                                        writeSPHDF5(FID,1,1,1,tau_D,tau_F,U);
                                                                                                    end
                                                                                                    
                                                                                                    %%%%%%% data sampling
                                                                                                    sample_pop = 1;
                                                                                                    writePopStatsRecordHDF5(FID, sample_pop);
                                                                                                    writePopStatsRecordHDF5(FID, 2);
                                                                                                    for pop_ind_pre = 1:Num_pop
                                                                                                        pop_ind_post = sample_pop;
                                                                                                        if pop_ind_pre == Num_pop
                                                                                                            syn_type = 2;
                                                                                                        else
                                                                                                            syn_type = 1;
                                                                                                        end
                                                                                                        
                                                                                                        % writeSynSampling(FID, pop_ind_pre, pop_ind_post, syn_type, sample_neurons, sample_steps)
                                                                                                        writeSynStatsRecordHDF5(FID, pop_ind_pre, pop_ind_post, syn_type)
                                                                                                    end
                                                                                                    
                                                                                                    %%% record synapse sample, especially for u&x
                                                                                                    writeSynSampling(FID,1,1,1,sample_neurons,ones(1, step_tot))
                                                                                                    
                                                                                                    %%% record neuron sample,especially for V&I
                                                                                                    % sample_steps = zeros(1,step_tot);
                                                                                                    % sample_steps(2*sec:step_tot) = 1;
                                                                                                    % writeNeuronSamplingHDF5(FID, sample_pop, [0,0,1,1,0,0,0,0],1:N(1), sample_steps)
                                                                                                    % writeNeuronSamplingHDF5(FID, sample_pop, [1,0,1,1,0,0,1,0], sample_neurons, ones(1, step_tot) )
                                                                                                    % writeNeuronSamplingHDF5(FID, 2, [1,1,1,1,0,0,1,1], [1 100], ones(1, step_tot) )
                                                                                                    
                                                                                                    % Add LFP sampling
                                                                                                    [Lattice, ~] = lattice_nD(2, hw);
                                                                                                    LFP_neurons = [];
                                                                                                    %                                                                                             LFP_centre_x = linspace(-hw,hw,63);
                                                                                                    %                                                                                             LFP_centre_y = linspace(-hw,hw,63);
                                                                                                    %                                                                                             LFP_centre_x = LFP_centre_x([16 47 32 32]); % 4 items
                                                                                                    %                                                                                             LFP_centre_y = LFP_centre_y([16 16 47 32]);
                                                                                                    LFP_centre_x = linspace(-hw,hw,9); % E16:9  E100:21 E400:41 E1600:81
                                                                                                    LFP_centre_y = linspace(-hw,hw,9);
                                                                                                    LFP_centre_x = LFP_centre_x(2:2:8); % E16(2:2:8)  E100(2:2:20) E400(2:2:40) E1600(2:2:80)
                                                                                                    LFP_centre_y = LFP_centre_y(2:2:8);
                                                                                                    [LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
                                                                                                    LFP_centre_x = LFP_centre_x(:);
                                                                                                    LFP_centre_y = LFP_centre_y(:);
                                                                                                    for cc = 1:length(LFP_centre_x)
                                                                                                        dist = lattice_nD_find_dist(Lattice, hw, LFP_centre_x(cc) , LFP_centre_y(cc));
                                                                                                        gaus_tmp = 1/(LFP_range_sigma*sqrt(2*pi))*exp(-0.5*(dist/LFP_range_sigma).^2) .* double(dist <= LFP_range_sigma*2.5);
                                                                                                        LFP_neurons = [LFP_neurons; transpose(gaus_tmp(:))]; %#ok<AGROW>
                                                                                                    end
                                                                                                    writeLFPRecordHDF5(FID, 1, LFP_neurons);
                                                                                                    
                                                                                                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                                                                                                    % Record explanatory variables, also called "controlled variables"
                                                                                                    
                                                                                                    writeExplVarHDF5(FID, 'discard_transient', discard_transient, ...
                                                                                                        'loop_num', loop_num, ...
                                                                                                        'delay', delay, ...
                                                                                                        'g_mu', g_mu, ...
                                                                                                        'g_EI', g_EI, ...
                                                                                                        'g_IE', g_IE, ...
                                                                                                        'g_II', g_II, ...
                                                                                                        'g_ext', g_ext,...
                                                                                                        'SpikeFreqAapt', SpikeFreqAapt, ...
                                                                                                        'dg_K',dg_K, ...
                                                                                                        'tau_K',tau_K,...
                                                                                                        'rate_ext_E', rate_ext_E,...
                                                                                                        'rate_ext_I', rate_ext_I,...
                                                                                                        'P0_init', P0_init, ...
                                                                                                        'degree_CV', degree_CV,...
                                                                                                        'in_out_r', in_out_r, ...
                                                                                                        'tau_c_EE', tau_c_EE, ...
                                                                                                        'tau_c_IE', tau_c_IE, ...
                                                                                                        'tau_c_I', tau_c_I, ...
                                                                                                        'STD_on', STD_on, ...
                                                                                                        'SP_on', SP_on, ...
                                                                                                        'cn_scale_wire', cn_scale_wire, ...
                                                                                                        'cn_scale_weight', cn_scale_weight, ...
                                                                                                        'mu_p', mu_p,...
                                                                                                        's_p', s_p, ...
                                                                                                        'inh_STDP', inh_STDP, ...
                                                                                                        'deg_hybrid', deg_hybrid,...
                                                                                                        'LFP_range_sigma', LFP_range_sigma,...
                                                                                                        'g_balance',g_balance,...
                                                                                                        'decay_GABA',decay_GABA,...
                                                                                                        'decay_AMPA',decay_AMPA,...
                                                                                                        'tau_D',tau_D, ...
                                                                                                        'tau_F',tau_F, ...
                                                                                                        'U',U, ...
                                                                                                        'coefficient_K_weight',Kw,...
                                                                                                        'NumP',NumP,...
                                                                                                        'ring_wid',ring_wid,...
                                                                                                        'ring',ring);
                                                                                                    %                                                                                                 'backgroundcoe',backgroundcoe);
                                                                                                    
                                                                                                    %%% save in_degree and sample neuron data based on in_degree
                                                                                                    % save([sprintf('%04g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF'),...
                                                                                                    %     '_config_data.mat'], 'in_degree', 'out_degree', 'dist_IJ', 'iter_hist',...
                                                                                                    %     'K_ee_mean', 'LFP_centre_x', 'LFP_centre_y','IndC','StiNeu'); % ,'I_ee','J_ee','K_ee','D_ee');
                                                                                                    
%                                                                                                     save([sprintf('%04g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF'),...
%                                                                                                         '_config_data.mat'],'StiNeu');
                                                                                                    
                                                                                                    % append this file self into .ygin for future reference
                                                                                                    appendThisMatlabFileHDF5(FID)
                                                                                                    
                                                                                                    disp('Matlab pre-processing done.')
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end

% This function must be here!
function appendThisMatlabFileHDF5(FID)
% need to coopy and past the following code into the file to be appended!
text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,'/config/MATLAB/config.m',text,'WriteMode','append');
end