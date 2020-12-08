function CH4Fig1
figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'CH4Fig1', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

id = [3 1 1];
for i = 1:3
    dir_strut = dir('*_RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
    dir_strut2 = dir('*_config_data.mat');
    num_files2 = length(dir_strut2);
    files2 = cell(1,num_files2);
    for id_out = 1:num_files2
        files2{id_out} = dir_strut2(id_out).name;
    end
    R = load(files{id(i)});
    load(files2{1},'StiNeu')    
    switch i
        case 1
            subplot(2,2,i)
            text(-0.1,1,'A','Units', 'Normalized','FontSize',12)
            cd ../4ItemsSimultaneousV2/
        case 2
            subplot(2,2,i)
            text(-0.1,1,'B','Units', 'Normalized','FontSize',12)
            cd ../SameCircuits4ItemsSequentialV2/
        case 3
            subplot(2,2,[3 4])
            text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
    end
    RasterPlotYL2(R,StiNeu)
end

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig1 % this is the trick!!
end