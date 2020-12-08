function [Dfre,Duration,Start,End] = BurstDrifting % (R)
% Adjust from function TimeFrequencyCombination1.m
tic
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
fs = 1e3;
lowFreq = 30;
hiFreq = 80; % Hz
Wn = [lowFreq hiFreq]/(fs/2);
order = 4; % 4th order
[b,a] = butter(order/2,Wn,'bandpass');
no = num_files;
seg = 1;
dt = 1;
i = 1;
% dt = R.dt;
% step_tot = R.step_tot;
% seg_size = 2e4; % 2s
% seg_num = ceil(step_tot/seg_size);
% [no, ~] = size(R.LFP.LFP_gamma);
Dfre = cell(1,no);
Duration = cell(1,no);
Start = cell(1,no);
End = cell(1,no);
% No.3
% freqrange = [30 80];
% Fs = 1000/dt;
% fc = centfrq('cmor1.5-1');
% scalerange = fc./(freqrange*(1/Fs));
% scales = scalerange(end):0.5:scalerange(1);
% pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs);
% scales = R.LFP.wavelet.scales;
% pseudoFreq = R.LFP.wavelet.pseudoFreq; % pseudo-frequencies
% for i = 1:no
%     for seg = 1:seg_num
%         seg_ind = get_seg(step_tot, seg_size, seg);
fsTemporal = 1/(dt*1e-3); % sampling frequency (Hz)
%
%         x_tmp = R.LFP.LFP_gamma(i,seg_ind);
for id_out = [3:13:130]
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    FR = vec2mat(R.num_spikes{1},10);
    FR = sum(FR,2)'/3969*1e3;
    FRcertain = FR; % filter(b,a,FR);
    % No.1 & 2
    [wt,f,coi] = cwt(FRcertain,fs,'VoicesPerOctave',30);
    for j = 1:length(coi)
        ind = find(f<=coi(j));
        wt(ind,j) = NaN;
    end
    % No.2
    ind = find(f<30 | f >80);
    wt(ind,:) = [];
    f(ind) = [];
    % No.3
%     f = pseudoFreq;
%     coeffs_tmp = abs(cwt(FR,scales,'cmor1.5-1'))'; % high spatial resolution; comr1-4, high frequency resolution
%     CData = transpose(coeffs_tmp); % coeffs_tmp'
%     zscoreCoef = zscore(CData,1);
%     zscoreCoef = CData;
%     Y = prctile(zscoreCoef(:),95) ;
    %     binaryImage = zeros(size(zscoreCoef)) ;
    %     binaryImage(zscoreCoef>=Y) = 1 ;
    CData = abs(wt);
    Y = prctile(CData(:),95); % 95
    CData(CData < Y) = 0;
    GreyImage = CData;
    CData(CData>0) = 1;
    binaryImage = CData;

    CC = bwconncomp(binaryImage) ;
    %     S = regionprops(CC,'Centroid'); % Center of mass of the region
    %     centroids = cat(1, S.Centroid);
    
    %         cTime = centroids(validRegionIdx(:),1)/fsTemporal ;
    %     cFreq = pseudoFreq(round(centroids(validRegionIdx(:),2))) ;
    
    B = regionprops(CC,'BoundingBox'); % Smallest rectangle containing the region
    Boundary = cat(1, B.BoundingBox); % [ul_corner width], ul_corner specifies
    % the upper-left corner [x y z ...]. width specifies the width along each dimension [x_width y_width ...]
    stTime = (seg-1)*2 + Boundary(:,1)*dt*1e-3; % second
    duTime = Boundary(:,3)/fsTemporal ;
    endTime = stTime + duTime;
    cFreq = [];
    iRegion = 2;
    j = 0;
    while iRegion < length(stTime)
        if stTime(iRegion+1) > endTime(iRegion)  %(size(CC.PixelIdxList{iRegion},1)>200)
%             certain = GreyImage(CC.PixelIdxList{iRegion});
%             [I,~] = ind2sub(size(GreyImage),CC.PixelIdxList{iRegion}(find(certain == max(certain))));
            
            instantBinary = binaryImage(:,floor(stTime(iRegion)*1e3):ceil(endTime(iRegion)*1e3));
            instantPattern = GreyImage(:,floor(stTime(iRegion)*1e3):ceil(endTime(iRegion)*1e3));            
        else
            while iRegion+1+j <= length(stTime) && stTime(iRegion+1+j) <= endTime(iRegion+j)
                stTime(iRegion+1+j) = NaN;
                endTime(iRegion+j) = NaN;
                j = j + 1;
            end
            if iRegion+1+j > length(stTime)
                break
            end
%             IdxList = cellfun(@transpose,CC.PixelIdxList,'UniformOutput',false);
%             certain = GreyImage([IdxList{iRegion:iRegion+j}]);
%             for n = iRegion:iRegion+j
%                 [I,~] = ind2sub(size(GreyImage),CC.PixelIdxList{n}(find(GreyImage(CC.PixelIdxList{n}) == max(certain))));
%                 if ~isempty(I)
%                     break
%                 end
%             end
            instantBinary = binaryImage(:,floor(stTime(iRegion)*1e3):ceil(endTime(iRegion+j)*1e3));
            instantPattern = GreyImage(:,floor(stTime(iRegion)*1e3):ceil(endTime(iRegion+j)*1e3));
        end
        iRegion = iRegion + j + 1;
        j = 0;        
%         cFreq = [cFreq f(I(1))];
        
        S = regionprops(instantBinary,instantPattern,{'Centroid','WeightedCentroid'});
        centroids = cat(1, S.WeightedCentroid);
        cFreq = [cFreq f(round(centroids(2)))];
    end
    stTime = stTime(~isnan(stTime))';
    endTime = endTime(~isnan(endTime))';
    stTime = stTime(2:end-1);
    endTime = endTime(2:end-1);
    %         freqLower = min([length(pseudoFreq)*ones(length(validRegionIdx),1),...
    %             round( centroids(validRegionIdx(:),2)+boundary(validRegionIdx(:),4)/2)]') ;
    %         freqUpper = max([1*ones(length(validRegionIdx),1),...
    %             round( centroids(validRegionIdx(:),2)-boundary(validRegionIdx(:),4)/2)]') ;
    %         bwFreq = pseudoFreq(freqUpper ) - pseudoFreq(freqLower) ;
    Dfre{i} = [Dfre{i} cFreq];
    Duration{i} = [Duration{i} endTime-stTime]; % second
    Start{i} = [Start{i} stTime];
    End{i} = [End{i} endTime];
    i = i + 1;
    toc
end
%     end
% end
% save('PopulationFRBurstDriftFrequencyPeakNoOverlap.mat','Dfre','Duration','Start','End')
%%
% edgesF = 30:5:90;
% N_F = zeros(no,length(edgesF)-1);
no = 100;
Interval = cell(1,no);
for i = 1:no
    %     [N_F(i,:),~] = histcounts(Dfre{i},edgesF,'normalization','pdf');
    Interval{i} = Start{i}(2:end) - End{i}(1:end-1);
end
% f = mean(N_F);
% std_f = std(N_F);
%%
figure
subplot(1,3,1)
f = cell2mat(Dfre);
histogram(f,18,'normalization','Probability')
xlabel('Drift Frequency(Hz)')
ylabel('Probability')
subplot(1,3,2)
d = cell2mat(Duration);
histogram(d,'normalization','Probability')
xlabel('Duration(s)')
ylabel('Probability')
subplot(1,3,3)
interval = cell2mat(Interval);
histogram(interval,'normalization','Probability')
xlabel('Interval(s)')
ylabel('Probability')
% text(-0.1,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%%
mat = cell(1,3);
mat{1} = A;
mat{2} = B;
mat{3} = C;
for i = 1:3
data = mat{i}; % (C>0);
figure(i)
subplot(1,2,1)
[N,edge]=histcounts(data(:),50,'normalization','pdf');
f = fitdist(data','Exponential')
x = linspace(0,max(data),1000);
semilogy(x,pdf(f,x),'r--')
f = fitdist(data','Lognormal')
hold on
semilogy(x,pdf(f,x),'b--')
f = fitdist([data,-data]','Stable')
hold on
semilogy(x,2*pdf(f,x),'g--')
hold on
semilogy(edge(2:end),N,'k.')
legend('Exponential fit','Lognormal fit','Alpha Stable fit','Data')
switch i
    case 1
        xlabel('Duration(ms)')
    case 2
        xlabel('Interval(ms)')
    case 3
xlabel('Size') % Interval(ms)
end
ylabel('pdf')
text(10,0.01,sprintf('alpha = %0.2f',f.alpha))
subplot(1,2,2)
[N,edges] = histcounts(data,50,'normalization','pdf');
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'o')
switch i
    case 1
        xlabel('Duration(ms)')
    case 2
        xlabel('Interval(ms)')
    case 3
xlabel('Size') % Interval(ms)
end
ylabel('pdf')
% hold on;
% Y = Y(N > 0);
% N = N(N > 0);
% YY = Y(3:20);
% NN = N(3:20);
% v = polyfit(log10(YY),log10(NN),1);
% x = Y(2:end);
% y = 10^v(2)*x.^v(1);
% loglog(x,y,'LineWidth',1.5)
% str = ['\alpha = ',num2str(-v(1))];
% text(min(x),2*max(y),str)
end
end