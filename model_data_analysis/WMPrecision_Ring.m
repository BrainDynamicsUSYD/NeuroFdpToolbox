function WMPrecision_Ring(varargin)
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% load('StiNeu.mat')
bin = 40; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
cent=((hw*2+1)^2+1)/2;
ring=20;
ring_wid=8.2;
dist_cent=lattice_nD_find_dist(Lattice,hw,cent); 
[neuron_dist,IndP] = sort(dist_cent); 
ring_dim=[ring-(ring_wid/2) ring+(ring_wid/2)]';
for i=1:length(ring_dim)
    if i==1
        idx=neuron_dist<ring_dim(i);
    elseif i==2
        idx=neuron_dist>ring_dim(i);
    end
    IndP(idx)=0; 
end
RingNeu=nonzeros(IndP)';
loop_num = 0;
th = 40;
for id_out = 1:num_files
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAY_INDEX = varargin{1};
        if loop_num ~=  PBS_ARRAY_INDEX
            continue;
        end
    end
    
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    RecallA = [];
%     switch ceil(id_out/100)
%         case 1
%             StiNeu = StiNeu1;
%         case 2
%             StiNeu = StiNeu2;
%         case 3
%             StiNeu = StiNeu3;
%         case 4
%             StiNeu = StiNeu4;
%     end
%     for no = 1:length(StiNeu)
        r = sum(movsum(full(R.spike_hist{1}(RingNeu,:)),bin,2)); % StiNeu{no}
        candx = Lattice(RingNeu,1); % StiNeu{no}
        candy = Lattice(RingNeu,2);
        for i = 2.26e4:length(r)-bin/2 % 4.46e4 % length(r) %
            if r(i) >= th % 25
                % neurons within WM area
                spike_x_pos_o = repmat(candx,1,bin).*R.spike_hist{1}(RingNeu,i-bin/2+1:i+bin/2); % column vector
                spike_x_pos_o = spike_x_pos_o(R.spike_hist{1}(RingNeu,i-bin/2+1:i+bin/2)); % StiNeu{no}
                spike_y_pos_o = repmat(candy,1,bin).*R.spike_hist{1}(RingNeu,i-bin/2+1:i+bin/2);
                spike_y_pos_o = spike_y_pos_o(R.spike_hist{1}(RingNeu,i-bin/2+1:i+bin/2));
                an = atan2d(spike_y_pos_o,spike_x_pos_o)';
                RecallA = [RecallA an];
            end
        end
%     end
    save([sprintf('%04g-RecallAngleThreashold%d',id_out,th),'.mat'],'RecallA')
    disp('Done');
end
%% plot
dir_strut = dir('*RecallAngleThreashold40.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
RecallAMat = [];
for id_out = 301:400 %:num_files    
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    RecallAMat = [RecallAMat R.RecallA];
    histogram(RecallAMat)
end
histogram(RecallAMat)
xlabel('Angle')
ylabel('Count')
end