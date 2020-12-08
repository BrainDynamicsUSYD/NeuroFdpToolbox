function CFCbefore(varargin)
% Main Function for cross-frequency coupling analysis
% Adapted from the function find_MI_cfc.m
% Required sub-folders:
% Toolbox_CSC
% ToolNeuroPatt
% Data
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
% clear
% close all
% clc
% cd ..
% addpath(genpath([pwd,'/Toolbox_CSC']))
%
% load([pwd,'/Data/UtahArrayData/SimulationRawLFP2'])
tic;
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
fs = 1e4;

% Simulation 40s
% clear
% addpath(genpath([pwd,'/Toolbox_CSC']))
% Loop number for PBS array job
loop_num = 0;
for id_out = 1:num_files
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    load(files{id_out},'LFP');
    % load([pwd,'/0021-201902221239-34344_in_1550799692431_out_RYG'],'LFP')
    % sigIn = reshape(LFP.LFP{1}(:,1:10:end),4,4,[]) ;
    % fsTemporal = 1e3 ;
    rawLFP = LFP.LFP{1}(:,1:5e4);
    phaseBand = [1:0.5:14] ;
    ampBand = 30:5:80 ;
    phaseBandWid = 0.49 ;
    ampBandWid = 5 ;
    
    % Butterworth filter
    order = 4; % 4th order
    
    % calculate the modulation index to find the coupling globally
    optionMethod = 1 ;         % 1 for KL distance
    optionSur = 2 ;            % 2 for block surrogate
    
    MI_raw_average = zeros(length(phaseBand), length(ampBand)) ;
    MI_surr_average = zeros(length(phaseBand), length(ampBand)) ;
    
    for j = 1:size(rawLFP,1) % randperm(size(rawLFP,1),10)% 1:size(rawLFP,1)
        for i = 1:length(phaseBand)
            lowF = phaseBand(i) - phaseBandWid;
            higF = phaseBand(i) + phaseBandWid;
            Wn = [lowF higF]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            LFP_theta = filter(b,a,rawLFP(j,:)); %#ok<AGROW>
            % hilbert transformation & gaussian smoothing
            sigPhase(i,:) = angle(hilbert(LFP_theta)); %#ok<AGROW>
        end
        for i = 1:length(ampBand)
            lowF = ampBand(i) - ampBandWid;
            higF = ampBand(i) + ampBandWid;
            Wn = [lowF higF]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            LFP_gamma = filter(b,a,rawLFP(j,:)); %#ok<AGROW>
            % hilbert transformation & gaussian smoothing
            sigAmp(i,:) = abs(hilbert(LFP_gamma)); %#ok<AGROW>
        end
        [MI_raw, MI_surr, ~, ~] =  find_MI_cfc ...
            ( sigPhase, sigAmp, optionMethod, optionSur) ;
        MI_raw_average = MI_raw +  MI_raw_average;
        MI_surr_average = MI_surr + MI_surr_average;
    end
    
    
    %% Contour Plot
    
    plotMI = MI_surr_average/size(rawLFP,1);
    
    % subplot(2,2,4)
    contourf(phaseBand,ampBand,plotMI')
    % text(-0.18,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    xlabel('phase frequency (Hz)')
    ylabel('amplitude frequency (Hz)')
    title('Modulation Index')
    % caxis([0,4e-3])
    colorbar
    savefig(gcf,[sprintf('%04g', id_out+3),'CFCbefore.fig'])
end
toc;