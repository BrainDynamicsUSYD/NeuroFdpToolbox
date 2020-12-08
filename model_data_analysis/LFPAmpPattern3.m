function LFPAmpPattern3(varargin)
% Collect LFP amplitude pattern dynamics on 3D scale on 20*20 LFP
% grid/ (40*40)
% only consider one pattern
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 10;
%%% Loop number for PBS array job %%%
loop_num = 0;

for i = 1:num_files
    for P = [47]
        for minBurstTime = [30] % 15 ms
            % For PBS array job
            loop_num = loop_num + 1;
            if nargin ~= 0
                PBS_ARRAYID = varargin{1};
                if loop_num ~=  PBS_ARRAYID
                    continue;
                end
            end
            R = load(files{i},'LFP');
            s = size(R.LFP.LFP{1});
            dt = 0.1;
            fs = 1/(bin*dt*1e-3); % sampling frequency (Hz)
            LFP_grid = reshape(R.LFP.LFP{1}(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]);
            clear R
%             LFP_grid = LFP_grid - mean(LFP_grid,3); % demean
%             % Butterworth filter
%             order = 4; % 4th order
%             lowFreq = 30; % gamma band (default values for this function are 150-250 Hz)
%             hiFreq = 80;
%             Wn = [lowFreq  hiFreq]/(fs/2);
%             [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
%             LFP_grid = LFP_grid(:,:,1*fs+1:end); % discard transient dynamics
            %         gamma_power_grid = zeros(size(LFP_grid));
            for m = 1:sqrt(s(1))
                for j= 1:sqrt(s(1))
                    gamma_tmp = LFP_grid(m,j,:);
%                     gamma_tmp = filter(b,a,LFP_grid(m,j,:));
                    gamma_tmp = gamma_tmp(:);
%                     gamma_tmp = gamma_tmp(1*fs+1:end-1*fs); % cut the head & tail after filtering
                    hil_tmp = abs(hilbert( gamma_tmp));
                    gamma_amp_grid(m,j,:) = hil_tmp; % (0.2*fs+1:end-0.2*fs); % cut the head & tail after Hilbert transform
                end
            end
            clear LFP_grid
            
            % Collecting Process
            Y = prctile(gamma_amp_grid(:),P);
            s = size(gamma_amp_grid);
            sigBinary = zeros(s);
            for t = 1:s(3)
                A = gamma_amp_grid(:,:,t);
                GGrid = zeros(s(1:2));
                GGrid(A >= Y) = 1;
                sigBinary(:,:,t) = GGrid;
            end
            CC = bwconncomp(sigBinary,6);
            per = [1 1 0];
            CC = CC2periodic(CC,per);
            Area = regionprops(CC,'Area') ;
            count = 1;
            for iBurst = 1: size(CC.PixelIdxList,1)
                currentIdx = CC.PixelIdxList{iBurst} ;
                burstTimeEnd = floor((currentIdx(end)-1)/(s(1)*s(2))) +1 ;
                burstTimeStart = floor((currentIdx(1)-1)/(s(1)*s(2))) +1 ;
                Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;    %
                if Duration2(iBurst) < minBurstTime
                    continue
                end
                Duration(count) = burstTimeEnd-burstTimeStart+1 ;    % duration
                
                patternScale(count) = Area(iBurst).Area ;  % 4 total scale
                
                burstIdxTemp = zeros(s) ;
                burstIdxTemp(currentIdx) = 1 ;
                currentBurst = sigBinary.*gamma_amp_grid.*burstIdxTemp ;
                sumAmp(count) = sum(currentBurst(:)) ;     % sum of amplitude
                peakAmp(count) = max(currentBurst(:)) ;    % peak amplitude
                
                for iTime = burstTimeStart:burstTimeEnd % [burstTimeStart,burstTimeEnd]
                    if iTime == burstTimeStart
                        Centroids{count} = [];
                        WCentroids{count} = [];
                        instantScale{count} = [];
                        Width{count} = [];
                    end
                    instantPattern = currentBurst(:,:,iTime);
                    [~,I] = max(instantPattern(:));
                    [I_row, I_col] = ind2sub(s([1 2]),I);
                    instantPattern = circshift(instantPattern, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
                    instantBinary = burstIdxTemp(:,:,iTime);
                    instantBinary = circshift(instantBinary, [round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
                    instantScale{count} = [instantScale{count};iTime,sum(instantBinary(:))] ;  % instant scale
                    S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});                    
                    Centroids{count} = [Centroids{count};iTime cat(1, S.Centroid)-[round(s(1)/2)-I_row round(s(2)/2)-I_col]];
                    WCentroids{count} = [WCentroids{count};iTime cat(1, S.WeightedCentroid)-[round(s(1)/2)-I_row round(s(2)/2)-I_col]];
                    B = regionprops(instantBinary,'BoundingBox'); % Smallest rectangle containing the region
                    Boundary = cat(1, B.BoundingBox);
                    Width{count} = [Width{count} (Boundary(3)+Boundary(4))/2/40*600]; % um
                end
                if count>1
                    ind = find(WCentroids{count}(:,1) == burstTimeStart);
                    firstCentroidsLoc = WCentroids{count}(ind,2:3) ;
                    distCent(count) = Distance_xy(firstCentroidsLoc(1),firstCentroidsLoc(2),lastCentroidsLoc(1),lastCentroidsLoc(2),s(1));
                    
                    firstCentroidsTime = burstTimeStart ;
                    centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
                else
                    distCent = [];
                    centInterval = [];
                end
                ind = find(WCentroids{count}(:,1) == burstTimeEnd);
                lastCentroidsLoc = WCentroids{count}(ind,2:3);
                lastCentroidsTime = burstTimeEnd ;
                count = count+1 ;
            end
            saveFileName = ['3DBurstRawLFP',sprintf('%04g',i),...
                'minTime',num2str(minBurstTime),'SR',num2str(fs),'P',num2str(P),'.mat'] ;
            save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp',...
                'patternScale','instantScale','Duration','Duration2','distCent','centInterval','Width') ;
        end
    end
end
end