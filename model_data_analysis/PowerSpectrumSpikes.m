function PowerSpectrumSpikes
%% plot power spectrum based on firing rate
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
pop = 1;
n = 1e4;
f = 1000/n*(0:round(n/2)); % start from 0 Hz

%%%%%%%%%DON'T CHANGE%%%%%%%%%%
if pop == 1
    pop_type = 'E';
else
    pop_type = 'I';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 % :num_files
%     fprintf('Loading RYG.mat file %s...\n', files{i});
%     R = load(files{i});
    y = vec2mat(R.num_spikes{pop},10); % spikes in ms    
    y = sum(y,2);
    y = y'/3969*1e3; % Hz on single neuronR
    Y = fft(y,n);
    PYY(i,:) = Y.*conj(Y)/n;
    Pyy(i,:) = PYY(i,1:round(n/2) + 1); 
end
pyy = mean(Pyy,1);
err = std(Pyy);
window_size = 10; % 1 Hz
simple = tsmovavg(pyy,'s',window_size,2);
% window_size = 50; % 5 Hz
% simple = tsmovavg(simple,'s',window_size,2);
% err = tsmovavg(err,'s',window_size,2);
% err = tsmovavg(err,'s',window_size,2);
x = f(11:2000);
y = simple(11:2000);
% y1 = y-err(1:1000);
% y2 = y+err(1:1000);
% X = [x,fliplr(x)];
% Y= [y1,fliplr(y2)];
% ind = find(~isnan(Y));
% X = X(ind);
% Y = Y(ind);
% fill(X,Y,[0.9 0.9 0.9]);
% hold on;
plot(x,y)
text(-0.28,1.02,'C','Units', 'Normalized','FontSize',14,'FontWeight','bold')
xlabel('Frequency(Hz)')
ylabel('Power')

%% Power Spectral Density Estimates Using FFT

Fs = 1e3;
pop = 1;
% y = R.Analysis.Hz_t{1};
y = sum(vec2mat(R0.num_spikes{pop},10),2); % firing rate
y = y'/3969*1e3; % Hz on single neuron
N = length(y);
ydft = fft(y);
ydft = ydft(1:N/2+1);
psdy = (1/(Fs*N)) * abs(ydft).^2;
psdy(2:end-1) = 2*psdy(2:end-1);
freq = 0:Fs/length(y):Fs/2;
window_size = 10; % 1 Hz
psdy = tsmovavg(psdy,'s',window_size,2);
plot(freq(4:1000),10*log10(psdy(4:1000)))
% plot(log10(freq(50:round(end/2))),log10(psdy(50:round(end/2))))
% hold on
% plot(log10(freq(round(end/5):round(end/2))),log10(psdy(round(end/5):round(end/2))))
% xlim([0 200])
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
end
