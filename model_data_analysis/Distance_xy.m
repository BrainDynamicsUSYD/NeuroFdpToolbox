function d = Distance_xy(I,J,x,y,l)
% In a l*l grid,distance between [I,J] and [x,y] considering periodic boundary
% condition. (I,J could be vectors)
% It's simple and helpful.
d_x = I-x;
d_y = J-y;
d_x(abs(I-x) <= l/2) = d_x(abs(I-x) <= l/2);
d_x(abs(I-x) > l/2) = abs(d_x(abs(I-x) > l/2))-l;
d_y(abs(J-y) <= l/2) = d_y(abs(J-y) <= l/2);
d_y(abs(J-y) > l/2) = abs(d_y(abs(J-y) > l/2))-l;
d = sqrt(d_x.*d_x+d_y.*d_y);
end
