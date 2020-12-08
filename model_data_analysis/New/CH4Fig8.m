function CH4Fig8
figure_width = 11.4; % cm
figure_hight = 8; % cm
figure('NumberTitle','off','name', 'CH4Fig8', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

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
hw = 31;
fw = 2*hw+1;
[Lattice, ~] = lattice_nD(2, hw);
ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];
post_dist = lattice_nD_find_dist(Lattice,hw,1);

R = load(files{32});
load(files2{21},'StiNeu','IndC')
subplot(2,3,[1 2])
RasterPlotYL2(R,StiNeu)
text(-0.1,1,'A','Units', 'Normalized','FontSize',12)

subplot(2,3,3)
box on;
[r,~] = sort(post_dist);
r = r(R.ExplVar.local_population);
for k = 1:length(IndC)
    r_cos = Lattice(IndC(k),1) + r*cos(ang); % N400:11.2, N100:5.66 N200:8.1
    r_sin = Lattice(IndC(k),2) + r*sin(ang);
    plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'k', ...
        r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'k',...
        r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'k',...
        r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'k',...
        r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'k', ...
        r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'k', ...
        r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'k', ...
        r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'k', ...
        r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'k');
    hold on
end
xlim([-hw hw]);
ylim([-hw hw]);
xlabel('x','FontSize',10)
ylabel('y','FontSize',10)
set(gca,'xtick',[],'ytick',[]);

R = load(files{34});
load(files2{23},'StiNeu','IndC')
subplot(2,3,[4 5])
RasterPlotYL2(R,StiNeu)
text(-0.1,1,'B','Units', 'Normalized','FontSize',12)

subplot(2,3,6)
box on;
[r,~] = sort(post_dist);
r = r(R.ExplVar.local_population);
for k = 1:length(IndC)
    r_cos = Lattice(IndC(k),1) + r*cos(ang); % N400:11.2, N100:5.66 N200:8.1
    r_sin = Lattice(IndC(k),2) + r*sin(ang);
    plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'k', ...
        r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'k',...
        r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'k',...
        r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'k',...
        r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'k', ...
        r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'k', ...
        r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'k', ...
        r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'k', ...
        r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'k');
    hold on
end
xlim([-hw hw]);
ylim([-hw hw]);
xlabel('x','FontSize',10)
ylabel('y','FontSize',10)
set(gca,'xtick',[],'ytick',[]);

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH4Fig8 % this is the trick!!
end