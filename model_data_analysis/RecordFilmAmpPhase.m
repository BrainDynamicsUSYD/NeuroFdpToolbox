% function RecordFilmAmpPhase
% record film of instantaneous gamma band LFP amp+phase on 63*63
load('0003-201804201449-05982_in_1524200018506_0_neurosamp.mat','LFP_grid','gamma_phase_grid','gamma_power_grid');
fs = 1e4; % sampling frequency (Hz)
fw = 63;
% Butterworth filter
order = 4; % 4th order
lowFreq = 1; % gamma band (default values for this function are 150-250 Hz)
hiFreq = 3;
Wn = [lowFreq  hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass');
delta_amp_grid = zeros(size(LFP_grid));
delta_phase_grid = zeros(size(LFP_grid));
disp('Creating delta_power_grid...')
for i = 1:fw
    for j= 1:fw
        delta_tmp = filter(b,a,LFP_grid(i,j,:));
        delta_tmp = delta_tmp(:);
        delta_amp_grid(i,j,:) = abs(hilbert(delta_tmp));
        delta_phase_grid(i,j,:) = angle(hilbert(delta_tmp));
    end
end
vidObj = VideoWriter('LFPGammaDeltaAmpPhase.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
fig = figure;
gamma_power_grid = gamma_power_grid(:,:,1:10:end);
gamma_phase_grid = gamma_phase_grid(:,:,1:10:end);
delta_amp_grid = delta_amp_grid(:,:,1:10:end);
delta_phase_grid = delta_phase_grid(:,:,1:10:end);
climAmp = minmax(gamma_power_grid(:)');
climPha = minmax(gamma_phase_grid(:)');
climAmp2 = minmax(delta_amp_grid(:)');
climPha2 = minmax(delta_phase_grid(:)');
for t = 1: 1e3 % 20001
    subplot(2,2,1)
    imagesc(gamma_power_grid(:,:,t),climAmp)
    title('Gamma Amp')
    subplot(2,2,2)
    imagesc(gamma_phase_grid(:,:,t),climPha)
    ts = sprintf('Gamma Phase, t = %8.1f ms',t);
    title(ts);
    subplot(2,2,3)
    imagesc(delta_amp_grid(:,:,t),climAmp2)
    title('Delta Amp')
    subplot(2,2,4)
    imagesc(delta_phase_grid(:,:,t),climPha2)
    title('Delta Phase')
    F = getframe(fig);
    writeVideo(vidObj,F.cdata);
    pause(0.05);
end
close(gcf);
close(vidObj);
% end