% test the feature of software transfer entropy
dbstop if error
bin = 300;

Opt.taux = -1; % arbitrary negative value>= tauy
Opt.trperm = 4;
Opt.method = 'dr';
Opt.bias = 'qe';
Opt.nt = 1; % number of calculating trails
Opt.testMode = 1;

% % Spike trains testing
% x = zeros(1,bin);
% x(randperm(bin,47)) = 1;
% y = zeros(1,bin);
% y(randperm(bin,47)) = 1;
% y(22:end) = x(1:end-21);
% X = x';
% Y = y';

% non-binary signal testing
t = 0:0.05:3*pi;
x = sin(t) + 1.2;
y = 1.2*rand(1,189);
% y(77:end) = x(1:end-76);
figure 
plot(t,x,t,y)
X = repmat(x,5,1)';
Y = repmat(y,5,1)';
Xsh = repmat(x(randperm(189)),5,1);
Ysh = repmat(y(randperm(189)),5,1);

te = zeros(1,bin);
nte = zeros(1,bin);
ntesh1 = zeros(1,bin);
ntesh2 = zeros(1,bin);
for j = 1:189
    Opt.tauy = -j;
    [TE,NTE] = transferentropy(X,Y,Opt,'TE','NTE');
    te(abs(j)) = nanmean(TE);
    nte(abs(j)) = nanmean(NTE);
    [TEsh,NTEsh] = transferentropy(Xsh',Ysh',Opt,'TE','NTE');
    ntesh1(abs(j)) = nanmean(TEsh);
    ntesh2(abs(j)) = nanmean(NTEsh);
end
figure
plot([1:bin],te,[1:bin],nte,[1:bin],ntesh1,'b--',[1:bin],ntesh2,'r--')
legend('te','nte','ntesh')
xlabel('tau_{delay}')
ylabel('transfer entropy')