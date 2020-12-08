function LFPAmpPattern5(varargin)
% Collect LFP pattern dynamics on 3D scale on 40*40 LFP
% using single electrode find burst method(baseline + n*s.d.)
% only consider one pattern
dbstop if error
dir_strut = dir('*RYG.mat');  % RYG % neurosamp
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

%%% Loop number for PBS array job %%%
loop_num = 0;
for i = 1:num_files
    for bin = [10]
        
        % For PBS array job
        loop_num = loop_num + 1;
        if nargin ~= 0
            PBS_ARRAYID = varargin{1};
            if loop_num ~=  PBS_ARRAYID
                continue;
            end
        end
        minBurstTime = 0; % 300/bin;
        fs = 1e4/bin; % Hz
        fprintf('Loading RYG.mat file %s...\n', files{i});
        R = load(files{i});
%         load(files{i},'LFP_grid');
%         R.LFP.LFP{1} = reshape(LFP_grid,63^2,[]);
%         clear LFP_grid
%         [R] = GetBurst2(R);
        sigBinary = R.LFP.GammaBurstEvent.is_burst(:,1:bin:end);
        s = size(sigBinary);
        sigBinary = flip(reshape(sigBinary,sqrt(s(1)),sqrt(s(1)),[])); % Actually no need to flip, sth legacy.
%         LFP_gamma_hilbert_abs = flip(reshape(R.LFP.LFP_gamma_hilbert_abs(:,1:bin:end),sqrt(s(1)),sqrt(s(1)),[]));
        LFP_gamma_hilbert_abs = flip(reshape(R.LFP.LFP_gamma_hilbert_abs,sqrt(s(1)),sqrt(s(1)),[]));
        clear R
        s = size(sigBinary);
        CC = bwconncomp(sigBinary,6);
        per = [1 1 0];
        CC = CC2periodic(CC,per);
        Area = regionprops(CC,'Area') ;
        count = 1;
%         disp('test 1')
        for iBurst = 1: size(CC.PixelIdxList,1)
            currentIdx = CC.PixelIdxList{iBurst} ;
            burstTimeEnd = floor((currentIdx(end)-1)/(s(1)*s(2))) +1 ;
            burstTimeStart = floor((currentIdx(1)-1)/(s(1)*s(2))) +1 ;
            Duration2(iBurst) =  burstTimeEnd-burstTimeStart+1 ;
            if Duration2(iBurst) < minBurstTime
                continue
            end
            Duration(count) = burstTimeEnd-burstTimeStart+1 ;
            
            patternScale(count) = Area(iBurst).Area ;
            
            burstIdxTemp = zeros(s) ;
            burstIdxTemp(currentIdx) = 1 ;
            currentBurst = sigBinary.*LFP_gamma_hilbert_abs.*burstIdxTemp ;
            sumAmp(count) = sum(currentBurst(:)) ;
            peakAmp(count) = max(currentBurst(:)) ;
%             disp('test 2')
            for iTime = burstTimeStart:burstTimeEnd % [burstTimeStart,burstTimeEnd]
                if iTime == burstTimeStart
                    Centroids{count} = [];
                    WCentroids{count} = [];
                    instantScale{count} = [];
                    Width{count} = [];
                end
%                 disp('test 3')
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
%                 disp('test 4')
                ind = find(WCentroids{count}(:,1) == burstTimeStart);
                firstCentroidsLoc = WCentroids{count}(ind,2:3) ;
                distCent(count) = Distance_xy(firstCentroidsLoc(1),firstCentroidsLoc(2),lastCentroidsLoc(1),lastCentroidsLoc(2),s(1));
                
                firstCentroidsTime = burstTimeStart ;
                centInterval(count) = firstCentroidsTime - lastCentroidsTime ;
            end
%             disp('test 5')
            ind = find(WCentroids{count}(:,1) == burstTimeEnd);
            lastCentroidsLoc = WCentroids{count}(ind,2:3);
            lastCentroidsTime = burstTimeEnd ;
            count = count+1 ;
        end
        saveFileName = ['3DBurst3',sprintf('%04g',i),...
            'minTime',num2str(minBurstTime),'SR',num2str(fs),'.mat'] ;
        save(saveFileName, 'Centroids','WCentroids', 'sumAmp','peakAmp',...
            'patternScale','instantScale','Duration','Duration2','distCent','centInterval','Width') ;
    end
end

% %% primary visualization
% % mat = reshape(R.LFP.GammaBurstEvent.is_burst,40,40,[]);
% mat = reshape(LFP.LFP{1},40,40,[]);
% lim = minmax(mat(:)');
% for i = 3911:10:length(mat) %
% %     imagesc(mat(:,:,i))
% %     caxis([0 1])
%     imagesc(mat(:,:,i),lim)
%     colorbar
%     ts = sprintf('t = %8.1f ms',0.1*i);
%     title(ts);
%     pause(0.01)
% end
end