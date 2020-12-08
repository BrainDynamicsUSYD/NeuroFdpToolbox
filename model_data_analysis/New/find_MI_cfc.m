function [ MI_raw, MI_surr, meanBinAmp, MI_surr_detail] = find_MI_cfc...
    ( sigPhase, sigAmp, optionMethod,optionSur)

% [MI_matrix_raw, MI_matrix_surr, MeanAmp, MI_surr] =  find_MI_cfc ... 
% ( sigIn, phaseBand, ampBand, optionMethod,optionSur,badChannels) ;
%
% This implements finding the modulation index of cross frequency coupling
% 
% Input:
%       sigPhase/Amp  (Freq X Time) Hilbert signals 
%       optionMethod  1 for Tort, 2 for Canolty, 3 for Ozkurt, 4 for PLV  
%       optionSur     1 for random phase and 2 for shifted phase
%
% Output:
%       MI_raw        (Phase X Amp) MI 
%       MI_surr       (Phase X Amp) MI with surrogate
%       MeanAmp       only available to 'Tort'
%
% Xian Long, Mar 19, 2018 @usyd. Supervisor: Pulin Gong
% xian.long@sydney.edu.au 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default settings
if ~exist('optionMethod','var')
    optionMethod = 1 ; 
end

if ~exist('optionSur','var')
    optionSur = 2 ;   
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSur = 200 ;
numPhase = size(sigPhase,1) ;
numAmp = size(sigAmp,1) ;
numBin = 20 ;
MI_raw = zeros(numPhase, numAmp) ;
MI_surr = zeros(size(MI_raw)) ;
meanBinAmp = zeros(numBin,numPhase, numAmp) ;
MI_surr_detail = zeros(numSur,numPhase, numAmp) ;

switch optionMethod
    case 1
        methodName = 'tort';
    case 2
        methodName = 'canolty' ;
    case 3
        methodName = 'ozkurt' ;
    case 4
        methodName = 'PLV' ;
end

if optionMethod ~= 1
    sigAmp = zscore(sigAmp) ;
end

for phaseNum = 1:numPhase
    for ampNum = 1:numAmp        
        [MI_raw(phaseNum,ampNum),MI_surr(phaseNum,ampNum), meanBinAmp(:,phaseNum,ampNum),...
            MI_surr_detail(:,phaseNum,ampNum)]...
            = CalcMI(sigPhase(phaseNum,:),sigAmp(ampNum,:),methodName,optionSur) ;
    end
    disp(['Complete percentage ',int2str(phaseNum/numPhase*100),'%'])
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MI_matrix_raw, MI_matrix_surr, MeanAmp, MI_surr] = ...
        CalcMI(phaseIn,ampIn,approach,surrogates)
% function CalcMI calculates the modulation index using four methods.
%
% author: Robert Seymour - Aston Brain Centre. July 2017
% edited by Xian Long
%%
phase = phaseIn ;
amp = ampIn ;

% Variable to hold MI for all trials
        MI_matrix_raw = [];
        
switch approach
    case 'tort'
        nbin = 20 ;
        [MI,MeanAmp] = calc_MI_tort(phase,amp,nbin);
        
    case 'ozkurt'
        [MI] = calc_MI_ozkurt(phase,amp);
        MeanAmp = zeros(20,1) ;
        
    case 'canolty'
        [MI] = calc_MI_canolty(phase,amp);
        MeanAmp = zeros(20,1) ;
        
    case 'PLV'
        [MI] = calc_MI_PLV(phase,amp);
        MeanAmp = zeros(20,1) ;
end

% Add the MI value to all other all other values
            MI_matrix_raw = MI;
            MI_matrix_surr = [] ;
            
% if strcmp(surrogates, 'yes')
 if surrogates==1 % random phase
     
    % Variable to surrogate MI
    MI_surr = [];
    
    % For each surrogate (surrently hard-coded for 200, could be changed)...
    for surr = 1:200
        % Get 2 random trial numbers
        trial_num = randperm(length(phaseIn),2);
        
        % Extract phase and amp info using hilbert transform
        % for different trials & shuffle phase
        phase=phaseIn(randperm(length(phaseIn))); % getting the phase
        amp = ampIn;
        
        % Switch PAC approach based on user input
        
        switch approach
            
            case 'tort'
                nbin = 20 ;
                [MI] = calc_MI_tort(phase,amp,nbin);
                
            case 'ozkurt'
                [MI] = calc_MI_ozkurt(phase,amp);
                
            case 'canolty'
                [MI] = calc_MI_canolty(phase,amp);
                
            case 'PLV'
                [MI] = calc_MI_PLV(phase,amp);
        end
        
        % Add this value to all other all other values
        MI_surr(surr) = MI;
    end
    
    
    % Subtract the mean of the surrogaates from the actual PAC
    % value and add this to the surrogate matrix
    MI_surr_normalised = MI_matrix_raw-mean(MI_surr);
    MI_matrix_surr = MI_surr_normalised;
 
 elseif surrogates==2 % random shifted phase, edited from Canolty et al.
     
     % Variable to surrogate MI
    MI_surr = [];
    
    numsurrogate = 200 ;
    % skip = ceil(length(phase).*rand(numsurrogate,1)); 

    numpoints = length(phase) ;
    % minskip = 1000 ;
    minskip= 200 ;
    maxskip = numpoints - minskip ;
    skip=ceil(numpoints.*rand(numsurrogate*2,1));       
    skip((skip>maxskip))=[];
    skip((skip<minskip))=[];
    skip=skip(1:numsurrogate,1);

    % For each surrogate (surrently hard-coded for 200, could be changed)...
    for surr = 1:numsurrogate
        % Get 2 random trial numbers
        amp=[ampIn(skip(surr):end), ampIn(1:skip(surr)-1)];
        
        phase=phaseIn;
        
        % Switch PAC approach based on user input
        
        switch approach
            
            case 'tort'
                nbin = 20 ;
                [MI] = calc_MI_tort(phase,amp,nbin);
                
            case 'ozkurt'
                [MI] = calc_MI_ozkurt(phase,amp);
                
            case 'canolty'
                [MI] = calc_MI_canolty(phase,amp);
                
            case 'PLV'
                [MI] = calc_MI_PLV(phase,amp);
        end
        
        % Add this value to all other all other values
        MI_surr(surr) = MI;
    end
    
    
    % Subtract the mean of the surrogaates from the actual PAC
    % value and add this to the surrogate matrix
    MI_surr_normalised = MI_matrix_raw-mean(MI_surr);
    MI_matrix_surr = MI_surr_normalised;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAC SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MI,MeanAmp] = calc_MI_tort(Phase,Amp,nbin)

% Apply Tort et al (2010) approach)
%nbin=18; % % we are breaking 0-360o in 18 bins, ie, each bin has 20o
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each bin
% (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi+(j-1)*winsize;
end

% now we compute the mean amplitude in each phase:
MeanAmp=zeros(1,nbin);
for j=1:nbin
    I = find(Phase <  position(j)+winsize & Phase >=  position(j));
    MeanAmp(j)=mean(Amp(I));
end

% The center of each bin (for plotting purposes) is
% position+winsize/2

% % Plot the result to see if there's any amplitude modulation
% if strcmp(diag, 'yes')
%     bar(10:20:720,[MeanAmp,MeanAmp]/sum(MeanAmp),'phase_freq')
%     xlim([0 720])
%     set(gca,'xtick',0:360:720)
%     xlabel('Phase (Deg)')
%     ylabel('Amplitude')
% end

% Quantify the amount of amp modulation by means of a
% normalized entropy index (Tort et al PNAS 2008):

MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
end

function [MI] = calc_MI_ozkurt(Phase,Amp)
% Apply the algorithm from Ozkurt et al., (2011)
N = length(Amp);
z = Amp.*exp(1i*Phase); % Get complex valued signal
MI = (1./sqrt(N)) * abs(mean(z)) / sqrt(mean(Amp.*Amp)); % Normalise
end

function [MI] = calc_MI_PLV(Phase,Amp)
% Apply PLV algorith, from Cohen et al., (2008)
amp_phase = angle(hilbert(detrend(Amp))); % Phase of amplitude envelope
MI = abs(mean(exp(1i*(Phase-amp_phase))));
end

function [MI] = calc_MI_canolty(Phase,Amp)
% Apply MVL algorith, from Canolty et al., (2006)
z = Amp.*exp(1i*Phase); % Get complex valued signal
MI = abs(mean(z));
end


end
