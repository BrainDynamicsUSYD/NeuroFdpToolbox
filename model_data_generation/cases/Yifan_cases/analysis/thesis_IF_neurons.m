% For most journals the figures should be 39 mm, 84 mm, 129 mm, or 174 mm wide and not higher than 234 mm.
figure_width = 15; % cm
figure_hight =  14; %cm
tick_fontsize = 7;
text_fontsize = 8;
anno_fontsize = 14;
anno_shift2 = [-0.07 0.08];


SV = 0.01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle','off','name','SPW-Spike lag analysis' ,'color','w', ...
  'units', 'centimeters', 'position', [0, 0, figure_width, figure_hight], ...
    'PaperSize', [figure_width, figure_hight] ,'PaperPositionMode','auto'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

dt = 0.1; %ms

Cm = 0.25; % nF
V_lk = -70.0;   % Leak reversal, -70.0
g_lk = 0.0167; % muS

V_rev = 0; % for AMPA
tau_r = 1.0; %ms
tau_d = 5; %ms

tau_V = 15;

exp_step_a = exp(-dt/tau_r);
exp_step_g = exp(-dt/tau_d);
exp_step_V = exp(-dt/tau_V);
t_step = 10^3;

a = zeros(1, t_step);
g = zeros(1, t_step);
I_a = zeros(1, t_step);
I_b = zeros(1, t_step);
V_a = zeros(1, t_step);
V_b = zeros(1, t_step);
V_c = zeros(1, t_step);

spikes = zeros(1,t_step);
spikes([100 ]) = 1;
spikes([100:100:500 ]) = 1;

ak = 0.6;
gk = 0.2;
k_I = 0.002;
for i = 2:t_step
    if spikes(i) == 1
        a(i) = a(i-1) + ak*(1-a(i-1));
        a(i) = a(i)*exp_step_a;
    else
        a(i) = a(i-1)*exp_step_a;
    end
  
    g(i) = g(i-1) + gk*a(i);
    g(i) = g(i)*exp_step_g;
    
    I_a(i) = k_I*g(i)*80;
    
    
    I_b(i) = k_I*g(i)*(80-V_b(i-1));
    
    
    V_a(i) = V_a(i-1) + I_a(i);
    V_a(i) = V_a(i)*exp_step_V;
    
    V_b(i) = V_b(i-1) + I_b(i);
    V_b(i) = V_b(i)*exp_step_V;
    
    if spikes(i) == 1
        V_c(i) =  V_c(i-1) + 8;
        V_c(i) = V_c(i)*exp_step_V;
    else
        V_c(i) = V_c(i-1)*exp_step_V;
    end
end


t = (1:t_step)*dt/1000;
figure(1);
subaxis(7,1,7,'SV', SV)
hold on;
d = 2;
s = find(spikes) - d/dt;
plot([s' s']'*dt/1000, [zeros(size(s')), ones(size(s'))]','k')
xlim(minmax(t))

subaxis(7,1,[6],'SV', SV)
hold on
plot(t,a,'k')
axis off;
subaxis(7,1,[5],'SV', SV)
hold on
plot(t,g,'b')
axis off;
subaxis(7,1,[3 4 ],'SV', SV)
hold on
plot(t,I_a,'r')
plot(t,I_b,'b')
axis off;
subaxis(7,1,[1 2],'SV', SV)
hold on
plot(t,V_a,'r')
plot(t,V_b,'b')
plot(t,V_c,'k')


