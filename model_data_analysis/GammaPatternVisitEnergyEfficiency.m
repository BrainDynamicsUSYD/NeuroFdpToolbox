% function VEF = GammaPatternVisitEnergyEfficiency
tic;
dir_strut = dir('*RYG.mat');  % RYG % neurosamp
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
hw = 31;
LFP_centre_x = linspace(-hw,hw,21); % E16:9  E100:21 E400:41 E1600:81
LFP_centre_y = linspace(-hw,hw,21);
LFP_centre_x = LFP_centre_x(2:2:20); % E16(2:2:8)  E100(2:2:20) E400(2:2:40) E1600(2:2:80)
LFP_centre_y = LFP_centre_y(2:2:20);
[LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
LFP_centre_x = LFP_centre_x(:);
LFP_centre_y = LFP_centre_y(:);
r = 0.005;
mat = zeros(1,130);
for id_out = 1:130
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    % spike
    R.step_tot = length(R.num_spikes{1});
    R.dt = 0.1;
    R.N = [3969 1000];
    [spike_hist_compressed,~] = find(R.spike_hist{1});
    R.spike_hist_compressed{1} = spike_hist_compressed';
    R = get_grid_firing_centreYL(R);
    
%     % LFP
%     th = 20; % prctile(R.LFP.LFP_gamma_hilbert_abs(:,500:end),95,'all');
%     interval = 1e3; % ms
%     count = 0;
    
    Ave = [];
    
%     % LFP
%     for i = 500:interval:length(R.LFP.LFP_gamma_hilbert_abs)-interval
%         ind = i:i+interval-1;
%         LFP_gamma = R.LFP.LFP_gamma_hilbert_abs(:,ind);
%         LFP_gamma(LFP_gamma<th) = 0;
%         WCentroid = [];
%         for j = 1:interval
%             LFP_gamma_hilbert_abs = LFP_gamma(:,j);
%             s = size(LFP_gamma_hilbert_abs);
%             LFP_gamma_hilbert_abs = reshape(LFP_gamma_hilbert_abs,sqrt(s(1)),sqrt(s(1))); % Actually no need to flip, sth legacy.
%             s = size(LFP_gamma_hilbert_abs);
%             sigBinary = LFP_gamma_hilbert_abs;
%             sigBinary(sigBinary > 0) = 1;
%             CC = bwconncomp(sigBinary,6);
%             per = [1 1];
%             CC = CC2periodic(CC,per);
%             if ~isempty(CC.PixelIdxList)
%                 if length(CC.PixelIdxList) > 1
%                     count = count + 1;
%                     [~,I] = sort(cellfun('length',CC.PixelIdxList),'descend');
%                     CC.PixelIdxList{1} = CC.PixelIdxList{I};
%                 end
%                 currentIdx = CC.PixelIdxList{1};
%                 burstIdxTemp = zeros(s) ;
%                 burstIdxTemp(currentIdx) = 1 ;
%                 sigBinary = sigBinary.*burstIdxTemp ;
%                 LFP_gamma_hilbert_abs = sigBinary.*LFP_gamma_hilbert_abs;
%                 [~,I] = max(LFP_gamma_hilbert_abs(:));
%                 [I_row, I_col] = ind2sub(s,I);
%                 LFP_gamma_hilbert_abs = circshift(LFP_gamma_hilbert_abs,[round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
%                 sigBinary = circshift(sigBinary,[round(s(1)/2)-I_row  round(s(2)/2)-I_col]);
%                 S = regionprops(sigBinary,LFP_gamma_hilbert_abs,'WeightedCentroid');
%                 WCentroid = [WCentroid; S.WeightedCentroid-[round(s(1)/2)-I_row round(s(2)/2)-I_col]];
%             end
%         end
%         if ~isempty(WCentroid)
%             [~,I] = min(pdist2([LFP_centre_x LFP_centre_y]+32,WCentroid));
%             Visit = sum(ismember(1:length(LFP_centre_x),I))/length(LFP_centre_x);
%         else
%             Visit = 0;
%         end
%         VisitEnergyEfficiency = Visit/(sum(R.num_spikes{1}(ind(1)*10:ind(end)*10))+r*(2*hw+1)^2*interval/50); % unit time 50 ms
%         Ave = [Ave VisitEnergyEfficiency];
%     end
    
    % spike
    interval = 2e4; % 0.1 ms
    for i = 5e3:interval:R.step_tot-interval
        ind = i:i+interval-1;
        WCentroid = R.grid.quick.centre(:,round((ind(1)-20)/10):round((ind(end)-20)/10)-1)';
        ind = ~isnan(WCentroid(:,1));
        WCentroid = WCentroid(ind,:);
        [~,I] = min(pdist2([LFP_centre_x LFP_centre_y],WCentroid));
        Visit = sum(ismember(1:length(LFP_centre_x),I))/length(LFP_centre_x);
        VisitEnergyEfficiency = Visit/(sum(R.num_spikes{1}(ind))+r*(2*hw+1)^2*interval/500); % unit time 50 ms
        Ave = [Ave VisitEnergyEfficiency];
    end
    
    mat(id_out) = mean(Ave);
    %     fprintf('\t Multiple Patterns:%d\n',count);
    toc;
end
mat = reshape(mat,13,[]);
VEF = mean(mat,2);
STD = std(mat,0,2);
errorbar(0.7:0.05:1.3,VEF,STD,'o--')
% plot(0.7:0.05:1.3,VEF,'o-')
xlabel('IE ratio')
ylabel('Searching esfficiency')
% end