%% read variables
% load('0005-201803070954-33175_in_1520377144068_0_neurosamp.mat','ripple_power_grid','ripple_fit')
load('0006-201804201449-05982_in_1524200018506_0_neurosamp.mat','gamma_power_grid')
% 'ripple_fit','ripple_fit_goodness','peak','gamma_fit'
minmumA = min(min(min(gamma_power_grid)));
maxmumA = max(max(max(gamma_power_grid)));
% minmumP = min(min(min(theta_phase_grid)));
% maxmumP = max(max(max(theta_phase_grid)));
threshold = exp(-1/2)*maxmumA; % corresponding sigma
[A_row, A_col, steps] = size(gamma_power_grid);
[x_grid, y_grid] = meshgrid(1:A_row, 1:A_col);

%% show the LFP pattern
vidObj = VideoWriter('3DGammaLFPAmp.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30; % number of frames to display per second
open(vidObj);
h = figure;
% set(gcf,'Position',[675 320 1603 641])
for i = 1:steps % 3650 %
    %     I_row = peak(i,1);
    %     I_col = peak(i,2);
    %     A_c = circshift(A, [round(A_row/2)-I_row  round(A_col/2)-I_col]);
    %     subplot(2,1,1)
    %     plot(ripple_fit{i}, [x_grid(:), y_grid(:)], A_c(:))
    %     subplot(3,2,4)
    %     subplot(1,2,1)
    %     imagesc(gamma_power_grid(:,:,i),[minmumA maxmumA])
    %     colorbar;
    %     subplot(1,2,2)
    %     imagesc(theta_phase_grid(:,:,i),[minmumP maxmumP])
    %     colorbar;
    surf(gamma_power_grid(:,:,i));
    zlim([minmumA maxmumA])
    ts = sprintf('Time = %8.1f ms', i*0.1);
    title(ts);
    %     text(-0.33,1.02,'E','Units', 'Normalized','FontSize',14,'FontWeight','bold')
    
    pause(0.02)
    F = getframe(h);
    writeVideo(vidObj,F.cdata);
    %     writeVideo(vidObj, getframe(gca));
        if i > 0.3e4
            break;
        end
end
close(gcf);
close(vidObj);

% vidTitle = 'AmpGammaPhaTheta' ;
% vidObj = VideoWriter(vidTitle);
% vidObj.FrameRate = 20 ;
% open(vidObj);
% fig=figure ;
% set(gcf,'Position',[675 362 1651 599])
%
% for i = 1e4:10:steps
%     subplot(1,2,1)
%     %colorMapSpec = pmkmp_new;
%     % Plot signal grid
%     imagesc(gamma_power_grid(:,:,i),[minmumA maxmumA])
%     set(gca,'YDir','normal')
%     title(['Time =', int2str(i*0.1),'ms'])
%     %colormap(gca)
%     colorbar
%
%
%     subplot(1,2,2)
%     imagesc(theta_phase_grid(:,:,i),[minmumP maxmumP])
%     set(gca,'YDir','normal')
%     title(['Time =', int2str(i*0.1),'ms'])
%     colorbar
%
%     writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
% %     pause(0.02)
%     if i > 1.5e4
%         break;
%     end
% end
% close(vidObj);

%% plot a series of moving
figure
countPlot = 1;
for i = (1:50:500) + 2130
    subplot(2,5,countPlot)
    countPlot = countPlot + 1;
    imagesc(gamma_power_grid(:,:,i),[minmumA maxmumA])
    %     title(['Time at ', int2str(i*0.1),'ms'])
    % colorbar
    % pause
    % cla
    % print(['Image',int2str(i*2/1000*1000),'ms'],'-depsc')
end

%% show the LFP pattern width distribution
% load('0001-201803070954-34384_in_1520377146275_0_neurosamp.mat','ripple_fit'); % default model
% load('0003-201803070954-33113_in_1520377128562_0_neurosamp.mat','ripple_fit'); % global normalized model
% load('0004-201803070954-33166_in_1520377135358_0_neurosamp.mat','ripple_fit'); % local normalized model
width = [];
ripple_fit = gamma_fit;
for i = 1:length(ripple_fit)
    if threshold > ripple_fit{i}.b
        if ripple_fit{i}.h > threshold-ripple_fit{i}.b
            d = 2*ripple_fit{i}.sigma*sqrt(2*log(ripple_fit{i}.h/(threshold-ripple_fit{i}.b)));
            width = [width d];
        end
    else
        d = 63;
        width = [width d];
    end
end
width = 600/63*width;
subplot(3,2,6)
histogram(width)
text(-0.28,1.02,'F','Units', 'Normalized','FontSize',14,'FontWeight','bold')
xlabel('Width(um)')
ylabel('Count')
% title('LFP pattern width distribution')

%% show membrane potential V
% load('0006-201804201449-05982_in_1524200018506_0_neurosamp.mat','V');
% R = load('0006-201804201449-05982_in_1524200018506_out_RYG.mat');
% S = full(R.spike_hist{1}(:,2e4:end));
% S = flip(permute(reshape(S,[63 63 20001]),[2 1 3]));
% V = flip(permute(reshape(V,[63 63 20001]),[2 1 3]));
% minmumV = min(min(min(V)));
% maxmumV = max(max(max(V)));
vidObj = VideoWriter('MembraneSpikes.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20; % number of frames to display per second
open(vidObj);
h = figure;
set(gcf,'Position',[675 320 1603 641])
for i = 1:20001
    ax1 = subplot(1,2,1);
    imagesc(V(:,:,i),[minmumV maxmumV])
    colormap(ax1,'parula')
    colorbar;
    ax2 = subplot(1,2,2);
    imagesc(S(:,:,i))
    colormap(ax2,'gray')
    ts = sprintf('Time = %8.1f ms', i*0.1);
    title(ts);
    pause(0.02)
    F = getframe(h);
    writeVideo(vidObj,F.cdata);
    %     writeVideo(vidObj, getframe(gca));
    if i > 0.3e4
        break;
    end
end
close(gcf);
close(vidObj);