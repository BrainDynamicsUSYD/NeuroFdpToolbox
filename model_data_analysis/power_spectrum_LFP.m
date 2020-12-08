function power_spectrum_LFP(stdin)
% plot power spectrum based on LFP

if nargin == 0
    dir_strut = dir('*RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for i = 1:num_files
        files{i} = dir_strut(i).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end

for i = 1 % :num_files % [1 12 17 18]
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    disp('done.\n');
    
    dt = R.dt;
    step_tot = R.step_tot;
    seg_size = 1e4; % 10s
    seg_num = ceil(step_tot/seg_size);
    [R] = get_IEI(R);
    [no, ~] = size(R.LFP.LFP_gamma);
    
    freqrange = [5 100];
    Fs = 1000/dt;
    fc = centfrq('cmor1.5-1');
    scalerange = fc./(freqrange*(1/Fs));
    scales = scalerange(end):0.5:scalerange(1);
    pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies
    tic;
    for j = 1:no
        for seg = 1:(2*seg_num - 1)
            
            seg_ind = ((seg - 1)*seg_size/2 + 1):((seg + 1)*seg_size/2);
            x_tmp = R.LFP.LFP_gamma(j,seg_ind);
            
            coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
            
            p(seg,:) = nanmean(coeffs_tmp);
        end
        power(j,:) = sum(p)/2;
        toc;
    end
    power = sum(power)/no;
    plot(pseudoFreq,power);
    xlabel('Frequency(Hz)');
    saveas(gca,'Power_LFP.pdf');
end
end