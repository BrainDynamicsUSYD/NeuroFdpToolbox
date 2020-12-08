function [Ns] = FindNeurons(Lattice,hw,num,varargin)
% find the closest num neurons to the center
if length(varargin) == 1
    [dist] = lattice_nD_find_dist(Lattice,hw,varargin{1});
else
    [dist] = lattice_nD_find_dist(Lattice,hw,varargin{1},varargin{2});
end
[~,Ns] = sort(dist);
% Ns = Ns(2:(1+num));
Ns = Ns(1:num);
end