function phasev = FreemanPattern(R,fre,mode)
% generate amplitude/phase pattern based on Freeman reserach
%  AM: root mean sqaure(sqrt(mean(x1^2+...+xn^2)))
% Ref: Freeman, Walter J., and John M. Barrie. "Analysis of spatial patterns
%     of phase in neocortical gamma EEGs in rabbit." Journal of neurophysiology 84.3 (2000): 1266-1278.
tic;
% pick burst time
% R0 = load('0001-201801030039-48409_in_1514900643054_out_RYG.mat');
% R = GetBurst2(R);
[no,step_tot] = size(R.LFP.LFP_broad);
% Btot = [];
% for i = 1 % :no
%     tot = find(R0.LFP.GammaBurstEvent.is_burst(i,:));
%     Btot = [Btot tot];
% end
% Btot = sort(unique(Btot));

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

% % create segment index for time points
% SegInd = cell(1,step_tot);
% for t = 1:step_tot-seg_size+gap
%     if ceil(t/gap) <= seg_size/gap
%         SegInd{t} = 1:ceil(t/gap);
%     else
%         SegInd{t} = ceil(t/gap)-seg_size/gap+1:ceil(t/gap);
%     end
% end
% for t = step_tot-seg_size+gap+1:step_tot
%     SegInd{t} = ceil(t/gap)-seg_size/gap+1:seg_num;
% end
% BurstNum = sort(unique([SegInd{Btot}]));

% FFT
Fs = 1e4; % Hz 0.1ms
n = 2^nextpow2(seg_size);
f = 0:(Fs/n):(Fs/2-Fs/n); % Frequency
ind = find(round(f)==fre);
if isempty(ind)
    disp('Frequency not match!')
end

l = [-pi pi]; %  [200 1600]
phasev = zeros(1,seg_num);
% vidObj = VideoWriter('ConePattern44HzDynamics2.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 3;
% open(vidObj);

% j = BurstNum(1);
j = 1;
for i = 1:seg_num % BurstNum 152
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
            
            %             % spatial filtering
            %             zeroPadMatrix = zeros(32,32);  % without windowing
            %             zeroPadMatrix(12:21,12:21) = dist;
            %             woHighPass = zeros(32,32);  % without windowing
            %             woWindFFT = fft2(zeroPadMatrix);
            %             woHighPass(3:30,3:30) = woWindFFT(3:30,3:30);
            %             woBandPass = woHighPass;
            %             woBandPass(12:21,12:21) = zeros(10,10);
            %             woSig = ifft2(woBandPass);
            %             dist = woSig(12:21,12:21);
            
            % nonlinear regression
            [xData, yData, zData] = prepareSurfaceData( [1:sqrt(no)], [1:sqrt(no)], dist);
            % Set up fittype and options.
            %             ft = fittype( 'poly22' );
            %             disp('Start fitting... \n')
            ft = fittype( @(x_c,y_c,z_c,a, x,y) a*sqrt((x-x_c).^2+(y-y_c).^2)+z_c,...
                'independent', {'x', 'y'},...
                'dependent', 'z');
            % Fit model to data.
            [fitr,~] = fit( [xData, yData], zData, ft,'StartPoint', [5 5 0 pi]);
            %             disp('Finish fitting... \n')
            x0 = meshgrid(1:10);
            y0 = meshgrid(1:10)';
            C = fitr.a*sqrt((x0-fitr.x_c).^2+(y0-fitr.y_c).^2)+fitr.z_c;
            C = reshape(flip(C),[1,no]);
            Rwf = (P - C)/P;
            x = meshgrid(-10:20) ;
            y = meshgrid(-10:20)' ;
            z = fitr.a*sqrt((x-fitr.x_c).^2+(y-fitr.y_c).^2)+fitr.z_c;
            %             z = fitresult.p00 + x.*fitresult.p10 + y.*fitresult.p20+...
            %                 fitresult.p20*x.^2 + fitresult.p11.*x.*y + fitresult.p02*y.^2 ;
            % temp(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
            
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
%     subplot(3,2,1)
%     contourf(dist)
%     colorbar
%     caxis(l)
%     text(-0.4,1.02,'A','Units', 'Normalized','FontSize',14,'FontWeight','bold')
%     subplot(3,2,3)
%     % pattern cone
%     contourf(z) % the coordinates for contourf is from 1~31
%     hold on;
%     line(11+[0 10 10],11+[0 0 10],'Color','red','LineStyle','-')
%     hold on;
%     line(11+[0 0 10],11+[0 10 10],'Color','red','LineStyle','-')
%     set(gca,'XtickLabel',[0:10:20]);
%     set(gca,'YtickLabel',[-5:5:20]);
%     %     imagesc(dist)
%     %     surf(dist)
%     %     zlim(l)
%     colorbar
%     caxis(l)
%     text(-0.4,1.02,'B','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    
    %         % dot trajectory
    %         if i == 1 || j == 1
    %             xc = fitr.x_c;
    %             yc = fitr.y_c;
    %         end
    %         line([0 10 10],[0 0 10],'Color','red','LineStyle','-')
    %         hold on;
    %         line([0 0 10],[0 10 10],'Color','red','LineStyle','-')
    %         xlim([-10 20])
    %         ylim([-10 20])
    %         hold on;
    %         h1 = plot(fitr.x_c,fitr.y_c, 'r>', 'MarkerSize', 8);
    %         hold on;
    %         if (fitr.x_c>=-10 && fitr.x_c<=20) || (fitr.y_c>=-10 && fitr.y_c<=20) % Distance_xy(x_tmp,y_tmp,x0,y0,fw) < 18
    %             plot([xc,fitr.x_c],[yc,fitr.y_c],'g')
    %             j = 0;
    %         else
    %             j = 1;
    %         end
    
    
    %     ts = sprintf('Time segment from %8.1f ms to %8.1f ms', (i-1)*2,(i-1)*2 +128);
    %     title(ts);
    %     pause(0.2)
    
    
    %     writeVideo(vidObj, getframe(gca));
    %     if i > 300
    %         break;
    %     end
    
    %     delete(h1);
    %     if j == 1
    %         delete(findobj(gca,'Type','line','Color','g'));
    %     end
    %     xc = fitr.x_c;
    %     yc = fitr.y_c;
end
subplot(3,2,5)
histogram(phasev,2e4);
xlim([0 8])
xlabel('Phase Velocity(M/sec)')
ylabel('Count')
text(-0.31,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
% close(gcf);
% close(vidObj);
% toc;
end