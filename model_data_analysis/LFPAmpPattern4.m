function LFPAmpPattern4(varargin)
% Collect LFP amplitude pattern dynamics on 3D scale on 63*63 LFP grid
% only consider one pattern
dir_strut = dir('*RYG.mat'); % '*_0_neurosamp.mat' '0*-ModifiedLFPGammaGrid.mat' '0*LocalSpikesGrid.mat'
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%% Loop number for PBS array job %%%
loop_num = 0;
% k = [7];
for i = 1:num_files
    for bin = [10]
        for P = [95]
            for minBurstTime = [0]/bin % 30 ms
                % For PBS array job
                loop_num = loop_num + 1;
                if nargin ~= 0
                    PBS_ARRAYID = varargin{1};
                    if loop_num ~=  PBS_ARRAYID
                        continue;
                    end
                end
                power_grid = load(files{i});
                power_grid = power_grid.LFP{1};
                % Butterworth filter
                order = 4; % 4th order
                lowFreq = 100; % ripple band (default values for this function are 150-250 Hz)
                hiFreq = 250;
                Wn = [lowFreq hiFreq]/(1e4/2);
                [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
                for j = 1:size(power_grid,1)
                    power_grid(j,:) = filter(b,a,power_grid(j,:)); 
                    power_grid(j,:) = abs(hilbert(power_grid(j,:)));
                end
                power_grid = reshape(power_grid,63,63,[]);
                %                 load(files{i},'ripple_power_grid')
                %                 power_grid = ripple_power_grid;
                %                 load(files{i},'gamma_power_grid')
                power_grid = power_grid(:,:,1:bin:end);
                fs = 1e4/bin;
                % Collecting Process
                Y = prctile(power_grid(:),P);
                s = size(power_grid);
                sigBinary = zeros(s);
                for t = 1:s(3)
                    A = power_grid(:,:,t);
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
                    currentBurst = sigBinary.*power_grid.*burstIdxTemp ;
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
                        Width{count} = [Width{count} (Boundary(3)+Boundary(4))/2/63*600]; % um
                    end
                    if count>1
                        ind = find(WCentroids{count}(:,1) == burstTimeStart);
                        firstCentroidsLoc = WCentroids{count}(ind,2:3) ;
                        distCent(count) = Distance_xy(firstCentroidsLoc(1),firstCentroidsLoc(2),lastCentroidsLoc(1),lastCentroidsLoc(2),s(1));
                        
                        firstCentroidsTime = burstTimeStart ;
                        centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
                    end
                    ind = find(WCentroids{count}(:,1) == burstTimeEnd);
                    lastCentroidsLoc = WCentroids{count}(ind,2:3);
                    lastCentroidsTime = burstTimeEnd ;
                    count = count+1 ;
                end
                saveFileName = ['3DBurstLFP',sprintf('%04g',i+1),...
                    'minTime',num2str(minBurstTime*bin/10),'SR',num2str(fs),'P',num2str(P),'.mat'] ;
                save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp',...
                    'patternScale','instantScale','Duration','Duration2','distCent','centInterval','Width') ;
            end
        end
    end
end
end