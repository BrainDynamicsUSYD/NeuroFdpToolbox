function SSCorrelation
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for i = 1:num_files
    files{i} = dir_strut(i).name;
end
%% plot a color image on frequency axis(x,y): comodulogram
lf = 0.01:99;
hf = 1:100;
fs = 1e4;
mat = zeros(length(lf));
for i = 1 % :num_files
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    % Butterworth filter
    order = 4; % 4th order
    for m = 1:length(lf)
        lowFreq = lf(m);
        hiFreq = hf(m);
        Wn = [lowFreq hiFreq]/(fs/2);
        [b,a] = butter(order/2,Wn,'bandpass');
        S1 = filter(b,a,R.LFP.LFP{1}(1,:));
        S1 = abs(hilbert(S1));
        for n = 1:length(lf)            
            lowFreq = lf(n);
            hiFreq = hf(n);
            Wn = [lowFreq hiFreq]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass');
            S2 = filter(b,a,R.LFP.LFP{1}(2,:));
            S2 = abs(hilbert(S2));
            Corr = corrcoef(S1,S2);
            mat(m,n) = Corr(2);
        end
    end
    mat = flip(mat);
    imagesc(mat)
    colorbar
    xlabel('Frequency(Hz)')
    ylabel('Frequency(Hz)')
end
%% cohernece/correlation map
for i = 1 % :num_files
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    no = size(R.LFP.LFP{1},1);
    mat = zeros(sqrt(no));
    ref = randperm(no,1);
    S1 = R.LFP.LFP_gamma(ref,:);
    for j = 1:no
        S2 = R.LFP.LFP_gamma(j,:);
        Corr = corrcoef(S1,S2);
        mat(j) = Corr(2);
    end
    [row,col] = ind2sub(size(mat),ref);
    imagesc(mat)
    colorbar
    title('Gamma LFP Correlation Map')
    hold on;
    plot(col,row,'k*', 'MarkerSize', 10)
end
end