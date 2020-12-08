% function WMCapacity
figure_width = 15.6; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'WMCapacity', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!
%%
dir_strut = dir('*_RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
bin = 40; % 40; % 4ms
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
C = zeros(1,num_files);
for NumP = 1:7
    switch NumP
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
    R = load(files{1+100*(NumP-1)});
    for i = 1:NumP
        dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
        LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
    end
    for id_out = 1:100
        for no = 1:NumP
            fprintf('Processing output file No.%d out of %d...\n', id_out+100*(NumP-1), num_files);
            fprintf('\t File name: %s\n', files{id_out+100*(NumP-1)});
            R = load(files{id_out+100*(NumP-1)});
            r = sum(movsum(full(R.spike_hist{1}(LoalNeu{no},2.26e4:end)),bin,2));
            if max(r) >= floor(0.3*length(LoalNeu{no})) %
                C(id_out+100*(NumP-1)) = C(id_out+100*(NumP-1)) + 1; % max(r)/length(LoalNeu{no});
            end
        end
    end
end
C = vec2mat(C,100);
%
% subplot(1,2,2)
cap = mean(C,2); % mean(C,2); % max(C,[],2);
STD = std(C,0,2);
errorbar(1:7,cap,STD,'o-','LineWidth',1.5)
% errorbar(1:7,cap,STD,zeros(size(STD)),'o-','LineWidth',1.5)
% plot(1:7,cap,'o-','LineWidth',1.5)
hold on
plot(1:7,1:7,'--','LineWidth',1)
xlabel('Implementing Items','fontSize',10)
ylabel('Capacity','fontSize',10)
%%
set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc WMCapacity
% end