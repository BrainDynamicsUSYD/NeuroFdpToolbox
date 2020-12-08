function K_ee = ModifyLocalEE(I_ee,J_ee,K_ee,Kw,d) % ,pr)
% This function work closely with the main function, the aim is to
% modify several local population connection strengths.
pool = 631:(630+d*63); % 1260
Ind = find(ismember(I_ee,pool) & ismember(J_ee,pool));
% Ind = randsample(Ind,round(pr*length(Ind)));
K_ee(Ind) = Kw*K_ee(Ind);
disp('Finishing modifying local EE.')

% %% create gradient bumps in the local path
% pool = [950 957 966 980 999];
% LP = [4 8 16 32 64];
% hw = 31;
% [Lattice,~] = lattice_nD(2,hw);
% for i = 1:length(pool)
%     post_dist = lattice_nD_find_dist(Lattice,hw,pool(i));
%     [~,IndP] = sort(post_dist);
%     pool = IndP(1:LP(i));
%     Ind = find(ismember(I_ee,pool) & ismember(J_ee,pool));
%     Ind = randsample(Ind,length(pool));
%     K_ee(Ind) = Kw*K_ee(Ind);
% end
% disp('Finishing modifying local EE.')

% %% create feedforward connection among local populations
% hw = 31;
% M = 100;
% Pff = 0.1;
% Prc = 0.1;
% sites = randperm((2*hw+1)^2-1,4);
% sites = [(2*hw+1)^2 sites];
% [Lattice,~] = lattice_nD(2,hw);
% l = length(I_ee);
% for i = 1:length(sites)
%     post_dist = lattice_nD_find_dist(Lattice,hw,sites(i));
%     [~,IndP] = sort(post_dist);
%     R = rand(M);
%     R(1:101:end) = 1;
%     if i == 1
%         pre = IndP(1:M);
%     else
%         post = IndP(1:M);
%         [row,col] = find(rand(M)<Pff);
%         I_ee = [I_ee;pre(row)];
%         J_ee = [J_ee;post(col)];
%         K_ee = [K_ee;K_ee(randperm(l,length(row)))];
%         pre = post;
%     end
% end
end