function schematicLF
% This simulation illustrates a fast implementation of three dimensional
% Brownian motion, the output is the Euclidean distance between initial
% and final positions.
% To calculate the mean value of T runs, run the following code in the 
% Command window :
%
% >>T=100;
% >> for n=1:T
% >>Brownianmotion;
% >>close;
% >>D(n)=d;
% >>end
% >>figure; plot(D),title(' Distance ');
% >>dd=mean(D)
% (c) Youssef Khmou, Applied mathematics, may 2015.
%%
rng(48,'twister')   %66
N=200;
x=cumsum(100*randn(1,N));
y=cumsum(100*randn(1,N));
% z=cumsum(randn(1,N));
close all
h=figure;
box on
h1= subplot(1,2,1)
% view(3);
set(gca,'GridLineStyle','--')
hold on;
% plot3(x,y,z,'LineWidth',1.5)
plot(x,y,'LineWidth',1.5)
%  axis([min(x) max(x) min(y) max(y) min(z) max(z)]);
axis([min(x) max(x) min(y) max(y)]);
axis([-9000 1000 -6000 2000]);

xlabel('x');
ylabel('y');
% zlabel('z');
% initial radius
% r0=[x(1) y(1) z(1)];
r0=[x(1) y(1)];
% final radius
% rf=[x(end) y(end) z(end)];
rf=[x(end) y(end) ];
plot(r0(1),r0(2),'or','MarkerFaceColor','r');
plot(rf(1),rf(2),'ok','MarkerFaceColor','k');
% Line of Sight between initial and final state
%xx=linspace(r0(1),rf(1),10);
%yy=linspace(r0(2),rf(2),10);
%zz=linspace(r0(3),rf(3),10);
%plot3(xx,yy,zz,'k--','LineWidth',2);
grid on;
% Distance
d=sqrt(sum((rf-r0).^2));
% Information=strcat('Three dimensional Brownian Motion, d=',num2str(d),' units');
% title(Information ,'FontWeight','bold');
% view(-109,58);

%
% figure;
h2=subplot(1,2,2)
rng(23,'twister')   %66
% x2 = cumsum(random('Stable',1.4,0,1,0,size(x))) ;
% y2 = cumsum(random('Stable',1.4,0,1,0,size(x))) ;
r = random('Stable',1.4,0,1,0,size(x)) ;
dir = rand(size(x)) ;
x2 = cumsum(100*r.*cos(dir)) ;
y2 = cumsum(100*r.*sin(dir)) ;
hold on
% plot3(x,y,z,'LineWidth',1.5)
plot(x2,y2,'LineWidth',1.5)
%  axis([min(x) max(x) min(y) max(y) min(z) max(z)]);
axis([min(x2) max(x2) min(y2) max(y2)]);
axis([-9000 1000 -6000 2000]);

xlabel('x');
ylabel('y');
% zlabel('z');
% initial radius
% r0=[x(1) y(1) z(1)];
r0=[x2(1) y2(1)];
% final radius
% rf=[x(end) y(end) z(end)];
rf=[x2(end) y2(end) ];
% plot3(r0(1),r0(2),r0(3),'or','MarkerFaceColor','r');
plot(r0(1),r0(2),'or','MarkerFaceColor','r');

% plot3(rf(1),rf(2),rf(3),'ok','MarkerFaceColor','k');
plot(rf(1),rf(2),'ok','MarkerFaceColor','k');

% Line of Sight between initial and final state
%xx=linspace(r0(1),rf(1),10);
%yy=linspace(r0(2),rf(2),10);
%zz=linspace(r0(3),rf(3),10);
%plot3(xx,yy,zz,'k--','LineWidth',2);
grid on;
% Distance
d=sqrt(sum((rf-r0).^2));
% Information=strcat('Three dimensional Levy walk, d=',num2str(d),' units');
% title(Information ,'FontWeight','bold');

%%
stepBr = sqrt(diff(x).^2 + diff(y).^2) ;
[x_n,n] = histcounts(stepBr,10,'normalization','probability') ;
figure;loglog( (n(1:end-1)+n(2:end))/2 ,x_n)
xlim([50 max(n)])

stepLe = sqrt(diff(x2).^2 + diff(y2).^2) ;
[x_n,n] = histcounts(stepLe,20,'normalization','probability') ;
figure;loglog( (n(1:end-1)+n(2:end))/2 ,x_n)
xlim([50 1000])

