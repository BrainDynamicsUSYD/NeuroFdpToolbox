function ind = IndexLattice(Lattice,I,J)
for i = 1:length(I)
    x = find(ismember(Lattice(:,1),I(i)));
    y = find(ismember(Lattice(:,2),J(i)));
    ind(i) = x(ismember(x,y));
%     display(i)
end
end