function [ I,J ] = Lattice2Lattice_rand( L_pre, L_post, p0 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


N_pre = length(L_pre(:,1));
N_post = length(L_post(:,1));
I = [];
J = [];
for i = 1:N_pre
    [~, ind] = sort( rand(N_post,1), 'ascend' );
    chosen_j = ind(1:poissrnd(N_post*p0));
    I = [I; i*ones(size(chosen_j))]; %#ok<AGROW>
    J = [J; chosen_j]; %#ok<AGROW>
end

end

