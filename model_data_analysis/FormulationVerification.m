% Formulation Verification
% It works on sample mat data.

dbstop if error
dir_strut0 = dir('*0_neurosamp.mat');
dir_strut1 = dir('*1_neurosamp.mat');
num_files = length(dir_strut0);
files0 = cell(1,num_files);
files1 = cell(1,num_files);
for i = 1:num_files
    files0{i} = dir_strut0(i).name;
    files1{i} = dir_strut1(i).name;
end
w = 40; % Hz
decay_E = 5e-3; % s
decay_I = [6e-3:1e-3:9e-3];
rise = 1e-3;
delay = 2e-3; % uniform 0~4 ms
AE = 5.58/7.38;
AI = 1/7.38;
SE = 1./sqrt((1+w.^2*decay_E^2).*(1+w.^2*rise^2));
SI = 1./sqrt((1+w.^2*decay_I.^2).*(1+w.^2*rise^2));
PhiE = w*delay + atan(w*rise) + atan(w*decay_E);
PhiI = w*delay + atan(w*rise) + atan(w*decay_I);
plot(decay_I,SI.*cos(PhiI),decay_I,SI.*sin(PhiI))
% plot(w,SE.*cos(PhiE),w,SE.*sin(PhiE),w,SI.*cos(PhiI),w,SI.*sin(PhiI))
% legend('Ecos','Esin','Icos','Isin')
% for i = 1:num_files % [1 12 17 18]
%     E = load(files0{i}); % 8 trails
%     I = load(files1{i}); % 2 trails
%     disp('Importing mat done.\n');
%     totE = E.I_AMPA + E.I_GABA + E.I_ext + E.I_leak;
%     totI = I.I_AMPA + I.I_GABA + I.I_ext + I.I_leak;
%     totE = mean(totE(:,(1e4+1):end),2);
%     totI = mean(totI(:,(1e4+1):end),2);
%     EE = mean(mean(E.I_AMPA(:,(1e4+1):end),2))/mean(totE);
%     IE = mean(mean(E.I_GABA(:,(1e4+1):end),2))/mean(totE);
%     EI = mean(mean(E.I_AMPA(:,(1e4+1):end),2))/mean(totI);
%     II = mean(mean(E.I_GABA(:,(1e4+1):end),2))/mean(totI);
%     Coe(i) = IE*EI-EE*II;
%     real(i) = SE*EE*AE*cos(PhiE) + SI*II*AI*cos(PhiI);
% %     disp(SE*EE*cos(PhiE))
% %     disp(SI*II*cos(PhiI))
% %     disp(SE*EE*sin(PhiE))
% %     disp(SI*II*sin(PhiI))
%     imagine(i) = SE*EE*AE*sin(PhiE) + SI*II*AI*sin(PhiI);
% end
