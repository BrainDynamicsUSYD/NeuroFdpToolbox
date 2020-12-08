 % Plot signal grid
 timeSlot = 600 ;
 signalGrid = abs(wvcfs(:,:, timeSlot)) ;
 sigLims = [min(abs(wvcfs(:))), max(abs(wvcfs(:)))];
 colorMapSpec = parula;
 vfColor = [1 1 1];
    imagesc(signalGrid, sigLims)
    colormap(gca, colorMapSpec)
    
    % Plot velocity field
    hold on
    vfScale = 1.0 ;
    vf = vfs(:,:,timeSlot) * vfScale;
    quiver(linspace(1, size(signalGrid,2), size(vf,2)), ...
        linspace(1, size(signalGrid,1), size(vf,1)), ...
        real(vf), imag(vf), 0, 'Color', vfColor)
    set(gca,'YDir','reverse');