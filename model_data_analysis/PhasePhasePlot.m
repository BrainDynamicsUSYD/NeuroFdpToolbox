% function PhasePhasePlot
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
Gamma = [];
Theta = [];
Delta = [];
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    LFP_gamma = R.LFP.LFP_gamma;
    no = size(LFP_gamma,1);
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 4;
    hiFreq = 10;
    fs = 1e4;
    Wn = [lowFreq  hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass');
    LFP_theta = zeros(size(LFP_gamma));
    for j = 1:no
        LFP_theta(j,:) = filter(b,a,R.LFP.LFP{1}(j,:));
    end
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 0.5;
    hiFreq = 4;
    fs = 1e4;
    Wn = [lowFreq  hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass');
    LFP_delta = zeros(size(LFP_gamma));
    for j = 1:no
        LFP_delta(j,:) = filter(b,a,R.LFP.LFP{1}(j,:));
    end
    for j = 1:no
        LFP_gamma(j,:) = angle(hilbert(LFP_gamma(j,:)));
        LFP_theta(j,:) = angle(hilbert(LFP_theta(j,:)));
        LFP_delta(j,:) = angle(hilbert(LFP_delta(j,:)));
    end
    Gamma = [Gamma;LFP_gamma(:)];
    Theta = [Theta;LFP_theta(:)];
    Delta = [Delta;LFP_delta(:)];
end
clear LFP_gamma LFP_theta LFP_delta
%%
dat = [Delta,Gamma];
n = hist3(dat,[90 90]);
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
h = pcolor(xb,yb,n1);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
%     oldcmap = colormap(gray);
colormap(flipud(colormap) );
colorbar
xlabel('Delta Phase(rad)')
ylabel('Gamma Phase(rad)')
% end