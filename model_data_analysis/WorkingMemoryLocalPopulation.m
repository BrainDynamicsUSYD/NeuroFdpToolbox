function [K_ee,StiNeu] = WorkingMemoryLocalPopulation(IndC,hw,Lattice_E,I_ee,J_ee,K_ee,Kw,LP)
% function [K_ee,Indegree,Pop] = WorkingMemoryLocalPopulation(NumP,I_ee,J_ee,K_ee,Kw,LP,pr)
% This function must work closely with the main function, the aim is to
% modify several local population connection strengths, i.e. make them
% identical and stronger.

% NumP = length(IndC);
% pool = cell(1,NumP);
% Kb = prctile(K_ee,50);
% Kp = prctile(K_ee,95);
% K_ee(1:end) = Kb;
% for i = 1:NumP
%     post_dist = lattice_nD_find_dist(Lattice_E,hw,IndC(i));
%     [~,IndP] = sort(post_dist);
%     pool{i} = IndP(1:LP)';
%     inn = zeros(1,LP);
%     for j= 1:LP
%         inn(j) = sum(ismember(I_ee,pool{i}) & ismember(J_ee,pool{i}(j)));
%     end
%     for j= 1:LP
%         Ind = find(ismember(I_ee,pool{i}) & ismember(J_ee,pool{i}(j)));
%         Ind = randsample(Ind,min(inn));
%         K_ee(Ind) = Kp;
%     end
% end
% All = 1:3969;
% other = setdiff(All,cell2mat(pool));
% for i = 1:NumP
%     Ind = find(ismember(I_ee,other));
%     Ind = randsample(Ind,round(0.1*length(Ind)));
%     K_ee(Ind) = Kp;
% end
% Indegree = [min(inn) Kb Kp];
% disp('Finishing adding local population for working memory.')

%% no local area for certain populations storing items
% Pop = cell(1,NumP);
% pool = randsample(1:3969,NumP*round(pr*LP));
% Pop{1} = pool(1:end/NumP);
% Pop{2} = pool((1+end/NumP):end);
% Indegree = [pool(1) 0 I_ee(J_ee==pool(1))'];
% for i = 1:NumP
%     Ind = find(ismember(I_ee,Pop{i}) & ismember(J_ee,Pop{i}));
%     K_ee(Ind) = Kw*K_ee(Ind);
% end
% disp('Finishing adding local populations for working memory.')

%% local poplulations are fixed on LFP electrodes
% And not every connection in the local population is increased
% At present, we use NumP = 2
% LP = 200; % neuron number in each local population
% pr = 0.6;
NumP = length(IndC);
pool = cell(1,NumP);
StiNeu = cell(1,NumP);
% IndAll = [];
% K = prctile(K_ee,90);
Kmean = mean(K_ee);
for i = 1:NumP
    post_dist = lattice_nD_find_dist(Lattice_E,hw,IndC(i));
    [~,IndP] = sort(post_dist);
    pool{i} = IndP(1:LP)';
    Ind = find(ismember(I_ee,pool{i}) & ismember(J_ee,pool{i}));
%     WMneurons{i} = unique(J_ee(Ind));
    K_ee(Ind) = Kw*max([K_ee(Ind) Kmean*ones(length(Ind),1)],[],2);
%     if i == 2
        StiNeu{i} = unique(I_ee(Ind));
%     end
%     IndAll = [IndAll;Ind];
end
% All = 1:3969;
% other = setdiff(All,cell2mat(pool));
% 
% Ind = find(ismember(I_ee,other));
% Ind = randsample(Ind,round(0.1*length(Ind)));
% K_ee(Ind) = Kw*K_ee(Ind);

% Indegree = [pool{1}(1) sum(ismember(I_ee,pool{1}) & ismember(J_ee,pool{1}(1)));...
%     pool{1}(end) sum(ismember(I_ee,pool{1}) & ismember(J_ee,pool{1}(end)))];
disp('Finishing adding local population for working memory.')

%% fixed local population + unchanged whole K
% % i.e. swap small K in local population with large K outside the local
% % population
% if NumP == 2
%     IndC = [513 2497];
% else
%     IndC = [513 2497 1994];
% end
% aInd = [];
% for i = 1:NumP
%     post_dist = lattice_nD_find_dist(Lattice_E,hw,IndC(i));
%     [~,IndP] = sort(post_dist);
%     pool = IndP(1:LP);
%     Ind = find(ismember(I_ee,pool) & ismember(J_ee,pool));
%     aInd = [aInd;Ind];
% end
% N = length(aInd);
% [B,I] = sort(K_ee);
% B = B(round(end/2):end);
% I = I(round(end/2):end);
% aInd = randsample(aInd,round(pr*length(aInd)));
% for i = 1:length(aInd)
%     tempK = K_ee(aInd(i));
%     tempI = [];
%     if tempK < B(1)
%         tempI = randsample(length(I),1);
%         K_ee(aInd(i)) = B(tempI);
%         K_ee(I(tempI)) = tempK;
%     else
%         tempI = find(I == aInd(i));
%     end
%     B(tempI) = [];
%     I(tempI) = [];
% end
% disp('Finishing adding local population for working memory.')

%% fixed local population + change other area K
% if NumP == 2
%     IndC = [513 2497];
% else
%     IndC = [513 2497 1994];
% end
% aInd = [];
% for i = 1:NumP
%     post_dist = lattice_nD_find_dist(Lattice_E,hw,IndC(i));
%     [~,IndP] = sort(post_dist);
%     pool = IndP(1:LP);
%     Ind = find(ismember(I_ee,pool) & ismember(J_ee,pool));
%     aInd = [aInd;Ind];
% end
% aInd = randsample(aInd,round(pr*length(aInd)));
% Ind = setdiff(1:length(I_ee),aInd);
% for i = Ind
%     K_ee(i) = 0.6*K_ee(i);
% end
% disp('Finishing adding local population for working memory.')
end