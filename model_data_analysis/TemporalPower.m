% function TemporalPower
figure_width = 11.4; % cm
figure_hight = 7.4; % cm
figure('NumberTitle','off','name', 'CH4Fig2', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 250; % ms
% Gamma_powerFR = cell(1,num_files-69); % num_files
% Gamma_powerLFP = cell(1,num_files-69);
xhere = 2e3:2e3:8e3; % ms
for id_out = 79 % :num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    try
        FR = vec2mat(R.num_spikes{1},10);
    catch
        FR = vec2mat(full(sum(R.spike_hist{1})),10);
    end
    FR = sum(FR,2)'/3969*1e3;
    FR = smooth(FR,bin);
    % Butterworth filter
    order = 4; % 4th order
    lowFreq_br = 30; % broad band (1-1000 Hz)
    hiFreq_br = 80;
    fs = 1e3;
    Wn = [lowFreq_br hiFreq_br]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    FR_gamma = filter(b,a,FR);
    FR_power = abs(hilbert(FR_gamma)).^2;
    FR_power = smooth(FR_power,bin);
%     t = 2e2:bin/2:length(FR_power)*0.6; % 0.1*R.step_tot;% R.step_tot;
    Gamma_powerFR = FR_power; % {id_out-69}
%     plot(t,FR_power(t),'.-','MarkerSize',10)
    
        try
            LFP_gamma = R.LFP.LFP_gamma;
            [no,~] = size(LFP_gamma);
        catch
            LFP = R.LFP.LFP{1};
            [no,~] = size(LFP);
            % Butterworth filter
            order = 4; % 4th order
            lowFreq_br = 30; % broad band (1-1000 Hz)
            hiFreq_br = 80;
            fs = 1e4;
            Wn = [lowFreq_br hiFreq_br]/(fs/2);
            [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
            LFP_gamma = zeros(size(LFP));
            for i = 1:no
                LFP_gamma(i,:) = filter(b,a,LFP(i,:));
            end
        end
        LFP_gamma_power = zeros(size(LFP_gamma));
        for i = 1:no
            LFP_gamma_power(i,:) = abs(hilbert(LFP_gamma(i,:))).^2;
        end
        Gamma_powerLFP = smooth(mean(LFP_gamma_power,1),bin*10); % {id_out-69}
%         t = 2e3:bin*5:length(Gamma_power{1})*0.6; % 0.6
%         plot(t*0.1,Gamma_power{id_out}(t),'.--','MarkerSize',10)
%         v = fit(t(3:80)'*0.1,Gamma_power{id_out}(t(3:80))','poly1');
%         hold on;
%         plot(t(3:80)*0.1,v.p1*t(3:80)*0.1+v.p2,'--')
    
    %     plot([xhere;xhere],[20*ones(1,4);160*ones(1,4)],'r-.')
    %     xlabel('Time(ms)')
    %     ylabel('Gamma Power')
    %     next = input('\t Next figure?');
    %     close all
end
% Gamma_powerFR = mean(cell2mat(Gamma_powerFR),2);
% Gamma_powerLFP = mean(cell2mat(Gamma_powerLFP),2);
t = 3e2:bin/2:length(Gamma_powerFR)*0.6; % *0.6
yyaxis left
plot(t*1e-3,Gamma_powerFR(t),'.-','MarkerSize',10)
ylim([2 9]*1e-4) % ([1 5]*1e-4) % 
ylabel('FR Gamma Power','fontsize',10)
t = 3e3:bin*5:length(Gamma_powerLFP)*0.6; % 0.6
yyaxis right
plot(t*0.1*1e-3,Gamma_powerLFP(t),'.-','MarkerSize',10)
hold on;
x = [xhere xhere+250 xhere+250 xhere]*1e-3;
y = [20*ones(1,8) 120*ones(1,8)]; % [10*ones(1,8) 80*ones(1,8)]; % 
v= [x' y'];
f = [1:4:16;2:4:16;3:4:16;4:4:16];
patch('Faces',f,'Vertices',v,'FaceColor',0.8*ones(1,3),'EdgeColor','none','FaceAlpha',0.3)
% plot([xhere;xhere],[20*ones(1,4);160*ones(1,4)],'r-.')
xlabel('Time(s)','fontsize',10)
ylabel('LFP Gamma Power','fontsize',10)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig2 % this is the trick!!
% end