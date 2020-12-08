function LFPAmpPattern2(varargin)
% Collect LFP amplitude pattern dynamics
% only consider one pattern
% adapt from LFPAmpPattern2.m and RecordAmpPattern2.m functions

dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 10;
Send = cell(1);
%%% Loop number for PBS array job %%%
loop_num = 0;

for i = 1:num_files
    for P = [95]
        for minBurstTime = [15] % 15 ms
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
            % Butterworth filter
            order = 4; % 4th order
            lowFreq = 5; % gamma band (default values for this function are 150-250 Hz)
            hiFreq = 100;
            Wn = [lowFreq  hiFreq]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); % The resulting bandpass and bandstop designs are of order 2n.
            gamma_amp_grid = zeros(size(LFP_grid));
            for m = 1:sqrt(s(1))
                for j= 1:sqrt(s(1))
                    gamma_tmp = filter(b,a,LFP_grid(m,j,:));
                    gamma_tmp = gamma_tmp(:);
                    hil_tmp = abs(hilbert( gamma_tmp));
                    gamma_amp_grid(m,j,:) = hil_tmp; % cut the head & tail after Hilbert transform
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
                ts(count) = burstTimeStart; % ms
                burstIdxTemp = zeros(s) ;
                burstIdxTemp(currentIdx) = 1 ;
                send = [];
                for iTime = burstTimeStart:burstTimeEnd
                    send = [send find(burstIdxTemp(:,:,iTime))'];
                end
                Send{count} = histcounts(send,0.5:(s(1)*s(2)+0.5));
                count = count+1 ;
            end
            saveFileName = ['3DBurst',sprintf('%04g',i),...
                'minTime',num2str(minBurstTime),'SR',num2str(fs),'P',num2str(P),'.mat'] ;
            save(saveFileName,'ts','Duration','Send') ;
        end
    end
end
% %% Data Analysis
% % pattern size and interval
% PS = patternScale(1:end-1);
% ESInterval = 0.1*centInterval(2:end);
% dat = [PS',ESInterval'];
% n = hist3(dat,[50 50]);
% n1 = n';
% n1(size(n,1) + 1, size(n,2) + 1) = 0;
% xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
% yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
% figure
% h = pcolor(xb,yb,n1);
% set(h, 'EdgeColor', 'none');
% h.ZData = ones(size(n1)) * -max(max(n));
% colormap(hot)
% oldcmap = colormap(gray);
% colormap( flipud(oldcmap) );
% colorbar
% [r,p] = corrcoef(PS',ESInterval');
% if p(2) > 10^(-4)
%     pvalue = ['p = ',num2str(p(2),'%.4f')];
% else
%     pvalue = ['p < 10^{-4}'];
% end
% tp = {['r = ',num2str(r(2),'%.4f')],pvalue};
% text(0.63,0.12,tp,'Units', 'Normalized','FontSize',8)
% xlabel('Pattern Size')
% ylabel('Interval among Bursts(ms)')
% figure
% [N,edges] = histcounts(patternScale,100);
% Y = (edges(1:end-1)+edges(2:end))/2;
% loglog(Y,N,'o')
% xlabel('Pattern Size')
% ylabel('Count')
end