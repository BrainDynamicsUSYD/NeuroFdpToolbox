figure_width = 11.4; % cm
figure_hight = 11.4; % cm
figure('NumberTitle','off','name', 'CH3Fig6', 'units', 'centimeters', ...
    'color','w', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight]); % this is the trick!

bin = 0.1; % ms
Size = [];
Duration = [];
vec = full(sum(R.spike_hist{1}));
mat = sum(reshape(vec,bin/R.dt,[]),1);
NonA = find(mat==0);
for i = 1:length(NonA)-1
    if NonA(i+1) - NonA(i) == 1
        continue
    end
    Size = [Size sum(mat(NonA(i):NonA(i+1)))];
    Duration = [Duration (NonA(i+1)-NonA(i)-1)*bin]; % ms
end
subplot(2,2,1)
edges = 10.^(linspace(log10(min(Size)),log10(max(Size)),6001));
N = histcounts(Size,edges,'Normalization', 'probability');
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'.')
hold on
Y = Y(N > 0);
N = N(N > 0);
range = 5:30;
v = polyfit(log10(Y(range)),log10(N(range)),1);
x = Y;
y = 10^v(2)*x.^v(1);
loglog(x,y,'linewidth',1.5)
str = ['\alpha = ',num2str(v(1))];
text(min(x),max(y),str)
xlabel('Size','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.2,1,'A','Units', 'Normalized','FontSize',12)

subplot(2,2,2) 
edges = 10.^(linspace(log10(min(Duration)),log10(max(Duration)),6001));
N = histcounts(Duration,edges,'Normalization', 'probability');
Y = (edges(1:end-1)+edges(2:end))/2;
loglog(Y,N,'.')
hold on
Y = Y(N > 0);
N = N(N > 0);
range = 1:15;
v = polyfit(log10(Y(range)),log10(N(range)),1);
x = Y;
y = 10^v(2)*x.^v(1);
loglog(x,y,'linewidth',1.5)
str = ['\alpha = ',num2str(v(1))];
text(min(x),max(y),str)
xlabel('Duration(ms)','fontsize',10)
ylabel('Probability','fontsize',10)
text(-0.2,1,'B','Units', 'Normalized','FontSize',12)

h = subplot(2,1,2);
p = get(h,'position');
Signal = R.grid.quick.jump_size;
[N,edges] = histcounts(Signal,100,'Normalization','pdf');
Y = (edges(1:end-1)+edges(2:end))/2;
pd = fitdist(Signal,'Stable')
y = pdf(pd,Y);

semilogx(Y,N,'.')
hold on
semilogx(Y,y,'r-.','LineWidth',1.5);
xlim([0 40])

xlabel('Step(unit grid)','fontsize',10)
ylabel('pdf','fontsize',10)
text(-0.1,1,'C','Units', 'Normalized','FontSize',12)
ax=axes('Position',[0.75 0.3 0.1 0.1],'Unit','normalize',...
    'parent',1);
loglog(ax,Y(Y>0),N(Y>0),'.')
hold on
loglog(ax,Y(Y>0),y(Y>0),'linewidth',1.5)

set(gcf, 'PaperPositionMode', 'auto'); % this is the trick!
print -depsc CH3Fig6 % this is the trick!!