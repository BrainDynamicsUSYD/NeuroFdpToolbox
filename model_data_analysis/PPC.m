function PPC %(varargin)
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% Loop number for PBS array job
loop_num = 0;
% PPCmat = cell(1,num_files);
for id_out = 1%:num_files
    % For PBS array job
    loop_num = loop_num + 1;
%     if nargin ~= 0
%         PBS_ARRAYID = varargin{1};
%         if loop_num ~=  PBS_ARRAYID
%             continue;
%         end
%     end
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    tic
    R = load(files{id_out});
    PPCvalue = 0; % zeros(1,51); % 30-80 Hz
    no = 1;
    for j = 30%:80
%         fprintf('Processing %d percent...\n', 2*(j-30));
        phase = cell(1,no);
        spikes = cell(1,no);
        Phase = cell(1,no);
        U = cell(1,no);
        P0 = 0;
        for i = 1:no
            HilbertOne = hilbert(R.LFP.LFP{1}(i,:))./abs(hilbert(R.LFP.LFP{1}(i,:)));
            [~,l1] = findpeaks(R.LFP.LFP{1}(i,:));
            tempICI = 0.1*(l1(2:end) - l1(1:end-1)); % ms
            freq = round(1e3./tempICI);
            ind = find(freq>=30 & freq<=80); % (freq == j);
            l2 = length(ind);
            for k = 1:l2
                phase{i} = [phase{i} HilbertOne(l1(ind(k)):l1(ind(k)+1))];
                spikes{i} = [spikes{i} R.num_spikes{1}(l1(ind(k)):l1(ind(k)+1))];
            end
            Phase{i} = repelem(phase{i},spikes{i});
            U{i} = [real(Phase{i});imag(Phase{i})];
        end
        N = length(Phase{i});
        for a = 1:N
            U1 = U{i}(:,a);
            for b = 1:N
                U2 = U{i}(:,b);
                P0 = P0 + dot(U1,U2);
            end
        end
        PPCvalue(j-29) = (P0-N)/(N*(N-1));
    end
    save([sprintf('%04g-', loop_num),'PPC_LFPNo1.mat'],'PPCvalue');
%     PPCmat{id_out} = PPCvalue;
    toc
end
% save('PPC.mat','PPCmat')
end