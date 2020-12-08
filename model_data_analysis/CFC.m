function CFC(StiNeu) % (varargin)
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
fs = 1e3;

% Simulation 40s
% clear
% addpath(genpath([pwd,'/Toolbox_CSC']))
% Loop number for PBS array job
% loop_num = 0;
for id_out = 22:num_files
    % For PBS array job
    %     loop_num = loop_num + 1;
    %     if nargin ~= 0
    %         PBS_ARRAYID = varargin{1};
    %         if loop_num ~=  PBS_ARRAYID
    %             continue;
    %         end
    %     end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    %     load(files{id_out},'LFP');
    % load([pwd,'/0021-201902221239-34344_in_1550799692431_out_RYG'],'LFP')
    % sigIn = reshape(LFP.LFP{1}(:,1:10:end),4,4,[]) ;
    % fsTemporal = 1e3 ;
    %     rawLFP = LFP.LFP{1}(1,:);
    R = load(files{id_out});
    Neuron = cat(1,StiNeu{:});
    num_spikes = sum(full(R.spike_hist{1}(Neuron,:)));
    FR = vec2mat(num_spikes,10);
    num = length(Neuron);
    FR = sum(FR,2)'/num*1e3;
    %     FR = vec2mat(R.num_spikes{1},10);
    %     FR = sum(FR,2)'/3969*1e3;
    rawLFP = FR(2.75e3+1:4.75e3); % (10e3+(1:10e3)); % FR(1:10e3); % FR(10e3+(1:10e3));
    rawLFP2 = FR(1:2e3);
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
    
    MI_raw_average2 = zeros(length(phaseBand), length(ampBand)) ;
    MI_surr_average2 = zeros(length(phaseBand), length(ampBand)) ;
    
    for j = 1:size(rawLFP,1) % randperm(size(rawLFP,1),10)% 1:size(rawLFP,1)
        for i = 1:length(phaseBand)
            lowF = phaseBand(i) - phaseBandWid;
            higF = phaseBand(i) + phaseBandWid;
            Wn = [lowF higF]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            LFP_theta = filter(b,a,rawLFP(j,:)); %#ok<AGROW>
            % hilbert transformation & gaussian smoothing
            sigPhase(i,:) = angle(hilbert(LFP_theta)); %#ok<AGROW>
            
            LFP_theta2 = filter(b,a,rawLFP2(j,:));
            sigPhase2(i,:) = angle(hilbert(LFP_theta2));
        end
        for i = 1:length(ampBand)
            lowF = ampBand(i) - ampBandWid;
            higF = ampBand(i) + ampBandWid;
            Wn = [lowF higF]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            LFP_gamma = filter(b,a,rawLFP(j,:)); %#ok<AGROW>
            % hilbert transformation & gaussian smoothing
            sigAmp(i,:) = abs(hilbert(LFP_gamma)); %#ok<AGROW>
            
            LFP_gamma2 = filter(b,a,rawLFP2(j,:)); %#ok<AGROW>
            sigAmp2(i,:) = abs(hilbert(LFP_gamma2));
        end
        [MI_raw, MI_surr, ~, ~] =  find_MI_cfc ...
            ( sigPhase, sigAmp, optionMethod, optionSur) ;
        MI_raw_average = MI_raw +  MI_raw_average;
        MI_surr_average = MI_surr + MI_surr_average;
        
        [MI_raw2, MI_surr2, ~, ~] =  find_MI_cfc ...
            ( sigPhase2, sigAmp2, optionMethod, optionSur) ;
        MI_raw_average2 = MI_raw2 +  MI_raw_average2;
        MI_surr_average2 = MI_surr2 + MI_surr_average2;
    end
    
    
    %% Contour Plot
    plotMI = MI_surr_average/size(rawLFP,1);
    
    subplot(1,2,2)
    contourf(phaseBand,ampBand,plotMI')
    % text(-0.18,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    xlabel('phase frequency (Hz)')
    ylabel('amplitude frequency (Hz)')
    title('Modulation Index')
    c2 = caxis;
    colorbar
    
    
    plotMI2 = MI_surr_average2/size(rawLFP2,1);
    
    subplot(1,2,1)
    contourf(phaseBand,ampBand,plotMI2')
    % text(-0.18,1.02,'D','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    xlabel('phase frequency (Hz)')
    ylabel('amplitude frequency (Hz)')
    title('Before-Modulation Index')
    c1 = caxis;
    c3 = [min([c1 c2]), max([c1 c2])];
    caxis(c3)
    colorbar
    
    subplot(1,2,2)
    caxis(c3)
    
    %     savefig(gcf,[sprintf('%04g', id_out+1),'PopFrCFC.fig'])
    
    %     figure
    %     bar(meanBinAmp(:,12,4))
    %     xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    %     ylabel('Phase Lock')
    %     savefig(gcf,[sprintf('%04g', id_out+1),'PopFrPhaseLock.fig'])
    next = input('\t Next figure?');
    close all
end
toc;