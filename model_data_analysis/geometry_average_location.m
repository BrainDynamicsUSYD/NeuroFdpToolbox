function [x,y] = geometry_average_location(x_combine,y_combine,hw)
% find the geometry average location of several points under spatial extended
% boundary condition
l = length(x_combine);
if l == 0
    x = [];
    y = [];
else
    x_standard = x_combine(1);
    y_standard = y_combine(1);
    x_combine = x_combine(2:l);
    y_combine = y_combine(2:l);
    x_combine(abs(x_combine - x_standard) <= hw) = x_combine(abs(x_combine - x_standard) <= hw);
    x_combine(abs(x_combine - x_standard) > hw) = x_combine(abs(x_combine - x_standard) > hw) + 2*hw*sign(x_standard - x_combine(abs(x_combine - x_standard) > hw));
    y_combine(abs(y_combine - y_standard) <= hw) = y_combine(abs(y_combine - y_standard) <= hw);
    y_combine(abs(y_combine - y_standard) > hw) = y_combine(abs(y_combine - y_standard) > hw) + 2*hw*sign(y_standard - y_combine(abs(y_combine - y_standard) > hw));
    x = (x_standard + sum(x_combine))/l;
    x(abs(x) <= hw) = x(abs(x) <= hw);
    x(abs(x) > hw) = x(abs(x) > hw) - 2*hw*sign(x(abs(x) > hw));
    y = (y_standard + sum(y_combine))/l;
    y(abs(y) <= hw) = y(abs(y) <= hw);
    y(abs(y) > hw) = y(abs(y) > hw) - 2*hw*sign(y(abs(y) > hw));
end
end