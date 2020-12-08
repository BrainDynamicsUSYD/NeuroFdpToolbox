% testing proper firing rate with tau_F
t_all = 1e5; % 10s
tau_D = 200; % ms
tau_F = 1500; % ms
dt = 0.1; % ms
U = 0.2;
N = 1;
u = 0.3*ones(N,t_all);
x = ones(N,t_all);
for t = 1:t_all-1
    nsp = rand(N,1)<2e-4;
    if t >= 5e4 && t<= 5.35e4
        nsp = rand(N,1)<0.25;
    end
    u(:,t+1) = u(:,t) + (U-u(:,t))/tau_F*dt + U*(1-u(:,t)).*nsp;
    x(:,t+1) = x(:,t) + (1-x(:,t))/tau_D*dt - u(:,t+1).*x(:,t).*nsp;
end
figure
subplot(2,1,1)
plot(0.1*(1:t_all),mean(u,1),0.1*(1:t_all),mean(x,1))
legend('u','x')
subplot(2,1,2)
plot(0.1*(1:t_all),mean(u,1).*mean(x,1))
% mean(u)

%% test gating variable dynamics
t_all = 1e5; % 10s
tau_d = 6; % ms
d = 4; % ms
dt = 0.1; % ms
N = 1;
s = zeros(N,t_all);
nspM = rand(N,t_all);
stepf = zeros(1,10); % last for 1ms
nspM(N,1:5e4-1) = nspM(N,1:5e4-1)<2e-4;
nspM(N,5e4:5.35e4) = nspM(N,5e4:5.35e4)<0.25;
nspM(N,5.35e4+1:end) = nspM(N,5.35e4+1:end)<2e-4;
for t = 1:t_all-1
    if t <= (d+1)/dt
        s(:,t+1) = s(:,t) - s(:,t)/tau_d*dt;
    else
        s(:,t+1) = s(:,t) - s(:,t)/tau_d*dt + dt*(1-s(:,t)).*logical(sum(nspM(:,t-(d+1)/dt:t-d/dt),2));
    end
end
figure
plot(0.1*(1:t_all),s)

%% test recover without spikes
t_all = 1e4; % 0.1ms
tau_D = 200; % ms
tau_F = 1500; % ms
dt = 0.1; % ms
u = zeros(1,t_all);
x = zeros(1,t_all);
U = 0.2;
u(1) = 0.92;
x(1) = 0.1;
for t = 1:t_all-1
    x(t+1) = (1-x(t))/tau_D*dt + x(t);
    u(t+1) = (U-u(t))/tau_F*dt + u(t);
end
plot(0.1*[1:t_all],x,'r.')
hold on;
plot(0.1*[1:t_all],u,'b.')
