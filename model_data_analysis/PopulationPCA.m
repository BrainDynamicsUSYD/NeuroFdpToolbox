% function PopulationPCA
% Ref: "Stable population coding for working memory coexists with
% heterogeneous neural dynamics in prefrontal cortex" PNAS
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
% dir_strut2 = dir('*_config_data.mat');
% load(dir_strut2.name,'StiNeu');
hw = 31;
[Lattice,~] = lattice_nD(2, hw);
NumP = 3;
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
R = load(files{46});
for i = 1:NumP
    dist = Distance_xy(Lattice(:,1),Lattice(:,2),Coor(1,i),Coor(2,i),2*hw+1); %calculates Euclidean distance between centre of lattice and node j in the lattice
    LoalNeu{i} = find(dist<=R.ExplVar.AreaR)';
end
N = min(cellfun('length',LoalNeu));
StiNeu = cell(1,NumP);
for i = 1:NumP
    StiNeu{i} = LoalNeu{i}(1:N);
end

bin = 2500; % 250ms bin
NumP = length(StiNeu);
r = cell(1,NumP);
period = 2.25e4+1:5.25e4;
% period = [5.25e4+1:5.6e4;10.25e4+1:10.6e4;15.25e4+1:15.6e4];
% period = [2.25e4+1:5.5e4;4.25e4+1:7.5e4;6.25e4+1:9.5e4;8.25e4+1:11.5e4];%...
%     10.25e4+1:10.59e4;12.25e4+1:12.59e4;14.25e4+1:14.59e4];
steps = length(period);
% steps = length(period)/bin;
Color = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 1 0;0 0 0];
for id_out = 1:NumP
%     fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
%     fprintf('\t File name: %s\n', files{id_out});
%     R = load(files{id_out});
    r{id_out} = movsum(full(R.spike_hist{1}(StiNeu{id_out},period)),bin,2)/bin*1e4; % Hz (id_out,:)
%     r{id_out} = squeeze(sum(reshape(full(R.spike_hist{1}(StiNeu{id_out},period)),length(StiNeu{id_out}),bin,[]),2))/bin*1e4; % Hz 
end
rmean = mean(cat(3,r{:}),3);   
z1 = zeros(NumP,steps);
z2 = zeros(NumP,steps);
X = zeros(NumP,length(StiNeu{1}));
for i = 1:NumP
    X(i,:) = mean(r{i},2) - mean(rmean,2);
end
[~,~,P] = svd(X);
% X = U*S*P', [U,S,P] = svd(X)
% X'*X = P*S'*S*P'
W1 = P(:,1);
W2 = P(:,2);
% vidObj = VideoWriter('0001-4PCA.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 6; % number of frames to display per second
% open(vidObj);
% fig = figure;
% xlabel('Stimulus PC1')
% ylabel('Stimulus PC2')
% zlabel('Time(ms)')
%
for t = bin:30:steps-bin % 1:steps % bin/2:steps-bin/2 % 
    %     for i = 1:NumP
    %         X(i,:) = r{i}(:,t) - rmean(:,t);
    %     end
    %     [~,~,P] = svd(X);
    %     % X = U*S*P', [U,S,P] = svd(X)
    %     % X'*X = P*S'*S*P'
    %     W1 = P(:,1);
    %     W2 = P(:,2);
    for i = 1:NumP
        z1(i,t) = W1'*(r{i}(:,t)-rmean(:,t));
        z2(i,t) = W2'*(r{i}(:,t)-rmean(:,t));
        plot3(z1(i,t),z2(i,t),0.1*t,'.','MarkerSize',4,'color',Color(i,:)) % t/steps*
%         hold on
%         plot(z1(i,t),z2(i,t),'.','MarkerSize',4,'color',Color(i,:))
        grid on
        drawnow
        hold on
        %         xlim([-200 200])
        %         ylim([-200 200])
%         F = getframe(fig);
%         writeVideo(vidObj,F.cdata);
    end
end
xlabel('Stimulus PC1')
ylabel('Stimulus PC2')
zlabel('Time(ms)')
% close(gcf);
% close(vidObj);
% end