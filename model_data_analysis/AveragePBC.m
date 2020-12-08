function mV = AveragePBC(Vec,l)
% average under periodic boundary condition, grid l*l
% Vec size n*2
mV = Vec/l*2*pi-pi;
mV = (angle(exp(1i*(mean(mV,1))))+pi)/(2*pi)*l;
end