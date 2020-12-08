function VVd = DistanceAmongMatrix(Mat1,Mat2)
% Mat [x1,y1;x2,y2;...]
% Vectors and distance [v1 v2 d;...]
l1 = size(Mat1,1);
l2 = size(Mat2,1);
VVd = zeros(l1*l2,5);
k = 1;
for i = 1:l1
    for j = 1:l2
        VVd(k,:) = [Mat1(i,:) Mat2(j,:) pdist([Mat1(i,:);Mat2(j,:)])];
        k = k + 1;
    end
end
% [dmin,I] = min(d);
% vec1 = Mat1(ceil(I/l2),:);
end