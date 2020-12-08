function FreemanPatternStable(varargin)
% generate amplitude/phase pattern based on Freeman reserach
% ajust from function FreemanPattern.m, calculate multiple frequencies in
% parallel
%  AM: root mean sqaure(sqrt(mean(x1^2+...+xn^2)))
% Ref: Freeman, Walter J., and John M. Barrie. "Analysis of spatial patterns
%     of phase in neocortical gamma EEGs in rabbit." Journal of neurophysiology 84.3 (2000): 1266-1278.
tic;
R = load('0001-201803041907-55109_in_1520151101934_out_RYG.mat');
[no,step_tot] = size(R.LFP.LFP_broad);
mode = 'phase';

% segments
seg_size = 1280; % 128 ms
gap = 20; % 2 ms
seg_num = floor((step_tot-seg_size)/gap) + 1;
Seg = zeros(seg_num,seg_size,no);
for j = 1:no
    for i = 1:seg_num
        Seg(i,:,j) = R.LFP.LFP_broad(j,1+20*(i-1):1280+20*(i-1));
    end
end

% FFT
Fs = 1e4; % Hz 0.1ms
n = 2^nextpow2(seg_size);
f = 0:(Fs/n):(Fs/2-Fs/n); % Frequency

% Loop number for PBS array job
loop_num = 0;
for fre = [29:5:59 63:5:83] % total 12
    ind = find(round(f)==fre);
    if isempty(ind)
        disp('Frequency not match!')
    end
    loop_num = loop_num + 1;
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    phasev = zeros(1,seg_num);
    Aone = [];
    Atwo = [];
    All = zeros(seg_num,2);
    for i = 1:seg_num % BurstNum
        seg = Seg(i,:,:); % 1*seg_size*no
        Y = fft(seg,n,2);
        switch mode
            case 'amplitude1'
                % FFT amplitude
                P2 = abs(Y/seg_size);
                P1 = P2(:,1:n/2+1,:);
                P1(:,2:end-1,:) = 2*P1(:,2:end-1,:);
                Amp = P1(:,1:n/2,:);
                dist = flip(reshape(Amp(:,ind,:),[sqrt(no) sqrt(no)]));
            case 'amplitude2'
                % Freeman AM
                seg = reshape(seg,[seg_size no]);
                RMS = sqrt(mean(sum(seg.^2),1));
                dist = flip(reshape(RMS,[sqrt(no) sqrt(no)]));
            case 'phase'
                P2 = angle(Y/seg_size);
                Phase = P2(:,1:n/2,:);
                dist = flip(reshape(Phase(:,ind,:),[sqrt(no) sqrt(no)]));
                P = reshape(Phase(:,ind,:),[1 no]);
                
                %                 % spatial filtering
                %                 zeroPadMatrix = zeros(32,32);  % without windowing
                %                 zeroPadMatrix(12:21,12:21) = dist;
                %                 woHighPass = zeros(32,32);  % without windowing
                %                 woWindFFT = fft2(zeroPadMatrix);
                %                 woHighPass(3:30,3:30) = woWindFFT(3:30,3:30);
                %                 woBandPass = woHighPass;
                %                 woBandPass(12:21,12:21) = zeros(10,10);
                %                 woSig = ifft2(woBandPass);
                %                 dist = woSig(12:21,12:21);
                
                % nonlinear regression
                [xData, yData, zData] = prepareSurfaceData( [1:sqrt(no)], [1:sqrt(no)], dist);
                ft = fittype( @(x_c,y_c,z_c,a, x,y) a*sqrt((x-x_c).^2+(y-y_c).^2)+z_c,...
                    'independent', {'x', 'y'},...
                    'dependent', 'z');
                
                % Fit model to data.
                [fitr,~] = fit( [xData, yData], zData, ft,'StartPoint', [5 5 0 pi]);
                x0 = meshgrid(1:10);
                y0 = meshgrid(1:10)';
                C = fitr.a*sqrt((x0-fitr.x_c).^2+(y0-fitr.y_c).^2)+fitr.z_c;
                C = reshape(flip(C),[1,no]);
                Rwf = (P - C)/P;
                if Rwf < 0.4
                    Aone = [Aone;fitr.x_c fitr.y_c Rwf -sign(fitr.a) i];
                    if Rwf < 0.2
                        Atwo = [Atwo;fitr.x_c fitr.y_c Rwf -sign(fitr.a) i];
                    end
                end
                All(i,:) = [fitr.x_c fitr.y_c];
                
                % calculate phase velocity
                dot1 = [5 5];
                p = polyfit([fitr.x_c fitr.y_c],dot1,1);
                if fitr.x_c < 5
                    dot2(1) = 6;
                elseif fitr.x_c > 5
                    dot2(1) = 4;
                else
                    dot1 = [6 6];
                    dot2(1) = 7;
                    p = polyfit([fitr.x_c fitr.y_c],dot1,1);
                end
                dot2(2) = polyval(p,dot2(1));
                phase1 = fitr.a*sqrt((dot1(1)-fitr.x_c).^2+(dot1(2)-fitr.y_c).^2)+fitr.z_c;
                phase2 = fitr.a*sqrt((dot2(1)-fitr.x_c).^2+(dot2(2)-fitr.y_c).^2)+fitr.z_c;
                slope = sqrt(sum((dot1-dot2).^2))*60e-3/abs(phase1-phase2); % mm/rad
                phasev(i) = slope*2*pi*f(ind)/1000;
        end
    end
    save([sprintf('0001-%04gNoFilter', fre),'Hz.mat'],'phasev','Aone','Atwo','All');
end
end