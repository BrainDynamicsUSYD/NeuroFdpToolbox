w = 1:200; % Hz
delay = 4e-3; % s
decay_E = 5e-3; % s
decay_I = 6e-3; % s
rise = 1e-3; % s
ratio = 0.7268;
S_E = 1./sqrt((1 + (w.^2).*decay_E^2).*(1 + (w.^2).*rise^2));
S_I = 1./sqrt((1 + w.^2*decay_I^2).*(1 + w.^2*rise^2));
Phi_E = w*delay + atan(w*rise) + atan(w*decay_E);
Phi_I = w*delay + atan(w*rise) + atan(w*decay_I);
y = sin(Phi_I) - S_E/S_I*ratio*sin(Phi_E);
plot(w,y)