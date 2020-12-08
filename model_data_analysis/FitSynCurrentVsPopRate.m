% fit the function v_alpha (t) = F{I_syn^alpha (t)}
% It seems that there is no specofic function I can work it out to fit it.
tic;
dir_strut1 = dir('*_0_neurosamp.mat'); % E
dir_strut2 = dir('*_1_neurosamp.mat'); % I
dir_strut3 = dir('*_out_RYG.mat'); 
num_files = length(dir_strut1);
files1 = cell(1,num_files);
files2 = cell(1,num_files);
files3 = cell(1,num_files);
Esyn = zeros(num_files,1e5);
Isyn = zeros(num_files,1e5);
vE = zeros(num_files,1e5);
vI = zeros(num_files,1e5);
j = 7;

for id_out = j % 1:num_files
    files1{id_out} = dir_strut1(id_out).name;
    files2{id_out} = dir_strut2(id_out).name;
    files3{id_out} = dir_strut3(id_out).name;
end
disp('Files loading done...')
for i = j % 1:num_files
    R = load(files1{i});
    Esyn(i,:) = nanmean(R.I_AMPA + R.I_GABA + R.I_ext);
end
disp('Esyn loading done...')
for i = j % 1:num_files
    R = load(files2{i});
    Isyn(i,:) = nanmean(R.I_AMPA + R.I_GABA + R.I_ext);
end
disp('Isyn loading done...')
for i = j % 1:num_files
    R = load(files3{i});
    vE(i,:) = R.Analysis.Hz_t{1};
    vI(i,:) = R.Analysis.Hz_t{2};
end
disp('vE/I loading done...')
h1 = plot(rand(100));
h2 = plot(rand(100));
xm1 = min(Esyn(j,5e4+1:1e5));
xm2 = max(Esyn(j,5e4+1:1e5));
ym1 = min(vE(j,5e4+1:1e5));
ym2 = max(vE(j,5e4+1:1e5));
xm3 = min(Isyn(j,5e4+1:1e5));
xm4 = max(Isyn(j,5e4+1:1e5));
ym3 = min(vI(j,5e4+1:1e5));
ym4 = max(vI(j,5e4+1:1e5));
for i = 5e4+1:1e5
    if i <= 5e4+100
        m = i-5e4;
    elseif mod(i,100) ~= 0
        m = mod(i,100);
        delete(h1(m));
        delete(h2(m));
    else
        m = 100;
        delete(h1(m));
        delete(h2(m));
    end
    subplot(1,2,1)
    h1(m) = plot(Esyn(j,i),vE(j,i),'bo');
    axis([xm1 xm2 ym1 ym2])
    drawnow    
    hold on
    subplot(1,2,2)
    h2(m) = plot(Isyn(j,i),vI(j,i),'r*');
    axis([xm3 xm4 ym3 ym4])
    drawnow
    hold on
    ts = sprintf('time = %8.1f ms', i*0.1);
    title(ts);
    pause(0.01)
end
toc;