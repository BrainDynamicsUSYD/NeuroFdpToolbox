% load('0001-201804201449-05982_in_1524200018506_0_neurosamp.mat', 'gamma_phase_grid');
% ft = fittype( @(x_c,y_c,z_c,a,b,c, x,y) ...
%     a*sqrt((angle(x./exp(1i*(x_c/63*2*pi-pi)))/(2*pi)*63).^2+(angle(y./exp(1i*(y_c/63*2*pi-pi)))/(2*pi)*63).^2)+z_c...
%     +2*pi*heaviside(sqrt((angle(x./exp(1i*(x_c/63*2*pi-pi)))/(2*pi)*63).^2+(angle(y./exp(1i*(y_c/63*2*pi-pi)))/(2*pi)*63).^2)-b)*(-1)^round(c),...
%     'independent', {'x', 'y'},...
%     'dependent', 'z');
vidObj = VideoWriter('PhasePattern.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 6; % number of frames to display per second
open(vidObj);
fig = figure;
set(gcf,'Position',[372 449 1213 420])
for i = 1:4e3 % length(gamma_phase_grid)
    %     [xData, yData, zData] = prepareSurfaceData( [1:63], [1:63],gamma_phase_grid(:,:,i));
    %     xData = exp(1i*(xData/63*2*pi-pi));
    %     yData = exp(1i*(yData/63*2*pi-pi));
    %     [fitr,~] = fit( [xData, yData], zData, ft,'StartPoint', [1 1 0 1 1 0]);
    %     x0 = meshgrid(1:63)/63*2*pi-pi;
    %     y0 = x0';
    %     z = fitr.a*sqrt((angle(x0./exp(1i*(fitr.x_c/63*2*pi-pi)))/(2*pi)*63).^2+(angle(y0./exp(1i*(fitr.y_c/63*2*pi-pi)))/(2*pi)*63).^2)...
    %         +fitr.z_c+2*pi*heaviside(sqrt((angle(x0./exp(1i*(fitr.x_c/63*2*pi-pi)))/(2*pi)*63).^2+...
    %         (angle(y0./exp(1i*(fitr.y_c/63*2*pi-pi)))/(2*pi)*63).^2)-fitr.b)*(-1)^round(fitr.c);
    subplot(1,2,1)
    phase_grid = gamma_phase_grid(:,:,i);
    phase_grid(phase_grid == -pi) = pi;
    imagesc(phase_grid,[-pi,pi])
    colorbar
    subplot(1,2,2)
    contourf(flipud(phase_grid))
    colorbar
    caxis([-pi pi])
    ts = sprintf('Time %8.1f ms', i*0.1);
    title(ts);
    F = getframe(fig);
    writeVideo(vidObj,F.cdata);
    pause(0.02)
end
close(gcf);
close(vidObj);