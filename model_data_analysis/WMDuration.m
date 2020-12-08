function Duration = WMDuration
% calculate Duration or Capacity
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
%%
dir_strut2 = dir('*_config_data.mat');
num_files2 = length(dir_strut2);
files2 = cell(1,num_files2);
for id_out = 1:num_files2
    files2{id_out} = dir_strut2(id_out).name;
end
bin = 40; % 4ms
%% Duration only item 1
Duration = zeros(1,num_files);
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
    load(files2{id_out},'StiNeu')
    r = sum(movsum(full(R.spike_hist{1}(StiNeu{1},:)),bin,2));
    Dind = find(r > 25);
    if Dind(end) > 2.25e4
        Duration(id_out) = Dind(end)*0.1-2.25e3;
    end
end
Duration = Duration*1e-3; % s
Duration = vec2mat(Duration,9);
dur = mean(Duration);
durSTD = std(Duration);
Fr = 0.6:0.1:1.4;
v2 = polyfit(log10(Fr(5:end)),log10(dur(5:end)),1);
x2 = Fr(4:end);
y2 = 10^v2(2)*x2.^v2(1);
errorbar(Fr,dur,durSTD,'o-.')
% set(gca,'YScale','log')
% hold on
% semilogy(x2,y2,'LineWidth',1.5)
xlabel('Scaled IE ratio')
ylabel('Working Memory Duration(s)')
% axes('Position',[.5 .5 .35 .35])
% box on
% errorbar(Fr,dur,durSTD,'o-.')
% hold on
% plot(x2,y2,'LineWidth',1.5)
% -v2(1)
% % errorbar(Fr,dur,durSTD,'o-')
ylim([0 25])
% xlabel('Scaled IE ratio')
% ylabel('Working Memory Duration(s)')
%% Average Duration
% Duration = cell(1,700);
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
    load(files2{id_out},'StiNeu')
    Duration{id_out+length(StiNeu)*100-100} = zeros(1,length(StiNeu));
    for i = 1:length(StiNeu) % length(StiNeu) % ceil(id_out/10)
        r = sum(movsum(full(R.spike_hist{1}(StiNeu{i},:)),bin,2));
        Dind = find(r > 25);
        if Dind(end) > 2.25e4
%             Duration{id_out+length(StiNeu)*100-100} = [Duration{id_out+length(StiNeu)*100-100} Dind(end)*0.1-2.25e3]; % +length(StiNeu)*100-100
            Duration{id_out+length(StiNeu)*100-100}(i) = Dind(end)*0.1-2.25e3;
        end
    end
end
%%
MeanDuration = cellfun(@mean,Duration)*1e-3; % s
% MeanDuration(isnan(MeanDuration)) = 0;
MeanDuration = vec2mat(MeanDuration,100);
dur = mean(MeanDuration,2);
STD = std(MeanDuration,0,2);
errorbar(1:7,dur,STD,'o-')
xlabel('Loading Items')
ylabel('Mean Recall Duration(s)')
%% Capacity
% Capacity = zeros(1,700);
inter = 1e4;
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    fprintf('\t File name: %s\n', files2{id_out}); % id_out-11
    load(files2{id_out},'StiNeu')    
    C = [];
    for t = 2.25e4+1:inter:(length(R.num_spikes{1})-inter)
        c = 0;        
        for i = 1:length(StiNeu)
            %         r = sum(movsum(full(R.spike_hist{1}(StiNeu{i},2.25e4+1:end)),bin,2));
            %         if max(r) > 25
            %             Capacity(id_out) = Capacity(id_out) + 1;
            %         end
            r = sum(movsum(full(R.spike_hist{1}(StiNeu{i},t:t+inter-1)),bin,2));
            if max(r) > 25
                c = c + 1;
            end
        end
        C = [C c];
    end
    Capacity(id_out+length(StiNeu)*100-100) = max(C);
end
%%
Capacity = vec2mat(Capacity,100);
cap = mean(Capacity,2);
STD = std(Capacity,0,2);
errorbar(1:7,cap,STD,'o-','LineWidth',1.5)
hold on
plot(1:7,1:7,'--','LineWidth',1)
xlabel('Implementing Items','fontSize',10)
ylabel('Capacity','fontSize',10)
%%
bin = 40; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
NumP = 3;
IE = 0.7:0.03:1.3;
switch 3 %NumP
    case 1
        Coor = [0;0];
    case 2
        Coor = [-16 15.5;-16 15.5];
    case 3
        Coor = [-10.5*sqrt(3) 10.5*sqrt(3) 0;-10.5 -10.5 21];
    case 4
        Coor = [-15.8 15.8 -15.8 15.8;-15.8 -15.8 15.8 15.8];
    case 5
        Coor = [-18.5 18.5 0 -18.5 18.5;-18.5 -18.5 0 18.5 18.5];
    case 6
        Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) -9.8*sqrt(3) 9.8*sqrt(3) 0;...
            -19.6 -9.8         -9.8        9.8          9.8         19.6];
    case 7
        Coor = [0 -9.8*sqrt(3) 9.8*sqrt(3) 0 -9.8*sqrt(3) 9.8*sqrt(3) 0;...
            -19.6 -9.8         -9.8        0 9.8          9.8         19.6];
end
LoalNeu = cell(1,NumP);
R = load(files{1});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
Duration = zeros(NumP,num_files);
for id_out = 1:num_files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = load(files{id_out});
    for i = 1:NumP
        r = sum(movsum(full(R.spike_hist{1}(LoalNeu{i},:)),bin,2));
        Dind = find(r > 25);
        if Dind(end) > 2.25e4
            Duration(i,id_out) = Dind(end)*0.1-2.25e3;
        end
    end
end
%%
duration = max(Duration)*1e-3; % s
duration = vec2mat(duration,length(IE));
dur = mean(duration);
durSTD = std(duration);
v2 = polyfit(log10(IE(9:end)),log10(dur(9:end)),1);
x2 = IE(7:end);
y2 = 10^v2(2)*x2.^v2(1);
errorbar(IE,dur,durSTD,'o-.')
% set(gca,'YScale','log')
hold on
semilogy(x2,y2,'LineWidth',1.5)
xlabel('Scaled IE ratio')
ylabel('Working Memory Duration(s)')
legend('WM duration','power-law fit')
end