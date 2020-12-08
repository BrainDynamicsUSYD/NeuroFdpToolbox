function RecordAmpPattern(varargin)
% 95% amplitude threshold for amplitude pattern
% cut head/tail version
close all;
clc;

%%% read all RYG.mat %%%
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

%%% Loop number for PBS array job %%%
loop_num = 0;

bin = 10;
P = 95;
ts = [];
for i = 1:num_files
    
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    fprintf('Loading RYG.mat file %s...', files{i});
    R = load(files{i});
    s = size(R.LFP.LFP{1});
    dt = 0.1;
    fs = 1/(bin*dt*1e-3); % sampling frequency (Hz)
    LFP_grid = flip(reshape(R.LFP.LFP{1}(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
    clear R
    LFP_grid = LFP_grid - mean(LFP_grid,3); % demean
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 5; % gamma band (default values for this function are 150-250 Hz)
    hiFreq = 100;
    Wn = [lowFreq  hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
    LFP_grid = LFP_grid(:,:,1*fs+1:end);
    amp_grid = zeros([sqrt(s(1)),sqrt(s(1)),s(2)/bin-3.4*fs]);
    for m = 1:sqrt(s(1))
        for j= 1:sqrt(s(1))
            tmp = filter(b,a,LFP_grid(m,j,:));
            tmp = tmp(:);
            tmp = tmp(1*fs+1:end-1*fs); % cut the head & tail after filtering
            hil_tmp = abs(hilbert(tmp));
            amp_grid(m,j,:) = hil_tmp(0.2*fs+1:end-0.2*fs); % cut the head & tail after Hilbert transform
        end
    end
    clear LFP_grid
    Y = prctile(amp_grid(:),P);
    s = size(amp_grid);
    for t = 1:s(3)
        A = amp_grid(:,:,t);
        GGrid = zeros(s(1:2));
        GGrid(A >= Y) = 1;
        CC = bwconncomp(GGrid);
        per = [1 1];
        CC = CC2periodic(CC,per);
        if isempty(CC.PixelIdxList)
            continue
        end
        ts = [ts t+2.2*fs];
%         [~,I] = max(A(:));
%         [I_row, I_col] = ind2sub(s([1 2]),I);
%         instantPattern = circshift(A, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
%         instantBinary = circshift(GGrid, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
%         S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});
    end 
    save(['AmpPatternLFPCut-',sprintf('%04g',loop_num),'.mat'],'ts');
end
end