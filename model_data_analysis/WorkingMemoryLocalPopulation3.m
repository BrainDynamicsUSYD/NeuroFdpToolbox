function [K_ee] = WorkingMemoryLocalPopulation3(hw,Lattice_E,I_ee,J_ee,K_ee,Coor,Kw,AreaR)
NumP = size(Coor,2);
Kmean = mean(K_ee);
for i = 1:NumP
    post_dist = Distance_xy(Lattice_E(:,1),Lattice_E(:,2),Coor(1,i),Coor(2,i),2*hw+1);
    pool = find(post_dist <= AreaR);
    Ind = find(ismember(I_ee,pool) & ismember(J_ee,pool));
    K_ee(Ind) = Kw*max([K_ee(Ind) Kmean*ones(length(Ind),1)],[],2);
end
disp('Finishing adding local population for working memory.')
end