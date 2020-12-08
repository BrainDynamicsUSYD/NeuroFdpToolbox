% A subplot format for distribution (2 variables)
[v,loop_number] = CollectVectorYG('grid','grid.jump_dist');
figure(1)
for i = 1:25
    subplot(5,5,i)
    histogram(v(loop_number >= (10*i - 9) & loop_number <= 10*i),100);
    axis([0 50 0 2000])
end
ax1 = axes('Position',[0 0 1 1],'Visible','off');% set nonvisible outside coordinate
axes(ax1);
title1 = 'Distribution of Jump Distance';
mark1 = 'Phi_I=0.8';
mark2 = 'Phi_I=0.9';
mark3 = 'Phi_I=1.0';
mark4 = 'Phi_I=1.1';
mark5 = 'Phi_I=1.2';
mark6 = 'Phi_E=0.8';
mark7 = 'Phi_E=0.9';
mark8 = 'Phi_E=1.0';
mark9 = 'Phi_E=1.1';
mark10 = 'Phi_E=1.2';
text(0.21,0.05,title1,'fontsize',14)
text(0.16,0.96,mark1)
text(0.32,0.96,mark2)
text(0.49,0.96,mark3)
text(0.65,0.96,mark4)
text(0.81,0.96,mark5)
text(0,0.87,mark6)
text(0,0.70,mark7)
text(0,0.52,mark8)
text(0,0.35,mark9)
text(0,0.17,mark10)
saveas(gcf,'Distribution_Jump_Distance.pdf');
figure(2)
for i = 1:25
    subplot(5,5,i)
    [N,X] = hist(v(loop_number >= (10*i - 9) & loop_number <= 10*i),100);
    semilogx(X,N);
    axis([10^-2 10^2 0 2000])
    %[N,X] = hist(v(loop_number > (20*i - 9) & loop_number < (20*i)));
    %loglog(X,N)
    %if N == zeros(1,10)
    %else
     %   f = fit(X',N','exp1');
      %  hold on
       % plot(f,X',N')
        %set(legend,'visible','off')
    %end
    %axis([10^0 10^2 10^(-1) 10^6])
end
ax2 = axes('Position',[0 0 1 1],'Visible','off');% set nonvisible outside coordinate
axes(ax2);
title2 = 'Semi-log Distribution of Jump Distance';
text(0.21,0.05,title2,'fontsize',14)
text(0.16,0.96,mark1)
text(0.32,0.96,mark2)
text(0.49,0.96,mark3)
text(0.65,0.96,mark4)
text(0.81,0.96,mark5)
text(0,0.87,mark6)
text(0,0.70,mark7)
text(0,0.52,mark8)
text(0,0.35,mark9)
text(0,0.17,mark10)
saveas(gcf,'Distribution_Jump_Distance_Semi-log.pdf');