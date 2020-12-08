function PostProcessYL2(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory
% adjust from PostProcesYL.m function
% Only read certain data to decrease the size of the data

% Prepare files
if nargin == 0
    dir_strut = dir('*out.h5');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end
% save figures
% save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for id_out = 1:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = ReadH5( files(id_out) ); % read .h5 file into matlab data struct
%     R = AnalyseYG(R); % do some simple analysis
%     fprintf('Starting getting Gamma...\n');
%     R = GetGamma(R{1});
%     R = GetBurst(R);
%     R = {R};
        R = R{1};
        
%         switch R.ExplVar.NumP
%             case 1
%                 Coor = [0;0];
%             case 2
%                 Coor = [-16 15.5;-16 15.5];
%             case 3
%                 Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
%             case 4
%                 Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
%             case 5
%                 Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
%             case 6
%                 Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%                     -19.6 -9.8         -9.8        9.8          9.8         19.6];
%             case 7
%                 Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
%                     -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
%         end
%         hw = 31;
%         [Lattice, ~] = lattice_nD(2, hw);
%         LoalNeu = cell(1,R.ExplVar.NumP);
%         for i = 1:R.ExplVar.NumP
%             dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
%             LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
%         end
%         RasterPlotYL2(R,LoalNeu)
%         saveas(gcf,[sprintf('%04g',63+R.ExplVar.NumP),'WMRasterPlot.eps'])
%         disp('Rasterplot Done');
        
        RL.stamp = R.stamp;
        RL.ExplVar = R.ExplVar;
%         RL.LFP = R.LFP.LFP;
%         RL.num_spikes = R.num_spikes;
        RL.spike_hist = R.spike_hist;
%         RL.reduced = R.reduced;
%         R = GetBurst2(R);
%         RL.LFP.GammaBurstEvent = R.LFP.GammaBurstEvent;
%         RL.LFP.LFP_gamma_hilbert_abs = R.LFP.LFP_gamma_hilbert_abs(:,1:10:end);
        clear R
        RL = {RL};
        SaveRYG(RL);
%     SaveRYG(R);
    disp('Done');
%     for ind = 1:length(R)
%         if isfield(R{ind}, 'samp_file')
%             Read_and_save_YGSamp(R{ind}.samp_file, R{ind});
%         end
%     end
%     RasterYL(R, save_fig); % generate raster plot for spiking history
end
end

