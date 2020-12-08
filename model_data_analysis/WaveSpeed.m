function WaveSpeed(varargin)
% calculate gamma wave speed
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end

% loop_num = 0;

velocity = [];
load('0001-201707201244-31484_config_data.mat')
% loop_num = loop_num + 1;

% For PBS array job
% if nargin ~= 0
%     PBS_ARRAYID = varargin{1};
%     if loop_num ~=  PBS_ARRAYID
%         continue;
%     end
% end
tic;
for i = 1:num_files % [1 12 17 18]
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    [Lattice, ~] = lattice_nD(2, 1.5);
    for i = 1:16
        d = lattice_nD_find_dist(Lattice,1.5,i);
        dist(i,:) = d';
        LFP_gamma_hilbert(i,:) = hilbert(R.LFP.LFP_gamma(i,:));
    end
    % dist = tril(dist);
    % Dis = unique(sort(dist));
    
    [send,receive] = find(dist);
    for i = 1e4:50:(R.step_tot-moments)
        
        fprintf('Calculating speed now...\n');
        
        phase1 = angle(LFP_gamma_hilbert(receive(ind),i));
        shift = angle(LFP_gamma_hilbert(send(ind),i)) - phase1;
        shift(shift < 0) = 2*pi + shift(shift < 0);
        shift(shift < 0) = pi + shift(shift < 0);
        k = i+1;
        phase2 = angle(LFP_gamma_hilbert(receive(ind),k));
        while sign(phase2) == sign(phase1)
            k = k + 1;
            phase2 = angle(LFP_gamma_hilbert(receive(ind),k));
        end
        p = phase2 - phase1;
        p(p<0) = p(p<0) + 3*pi;
        T = (k - i)/p*2*pi*R.dt; % ms
        speed = dist(send(ind),receive(ind))*150/(shift/2/pi*T); % um/ms = mm/s
        velocity = [velocity speed];
        
    end
end
histogram(velocity*1e-2,[1:100]) % cm/s
xlabel('Wave Speed(cm/s)')
ylabel('Number count')
saveas(gcf,['WaveSpeedDistribution',num2str(moments),'.pdf'])
toc;
end