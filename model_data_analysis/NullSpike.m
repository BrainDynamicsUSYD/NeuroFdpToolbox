function NullSpike(R)
GammaLFP = R.LFP.LFP_gamma;
GLFPPower = R.LFP.LFP_gamma_hilbert_abs.^2;
[no,tot] = size(GammaLFP);
t = 1e4+1:1e4+800; % 1s~4s
for i = 1:no
%     subplot(2,2,1)
%     plot(1e-4*t,GammaLFP(i,t))
%     xlabel('Time(s)')
%     ylabel('Gamma band LFP')
%     hold on;
%     subplot(2,2,2)
%     plot(1e-4*t,GLFPPower(i,t))
%     xlabel('Time(s)')
%     ylabel('Analytic Power')
%     hold on;
%     subplot(2,2,3)
%     plot(1e-4*t,log10(GLFPPower(i,t)));
%     xlabel('Time(s)')
%     ylabel('log_{10}Analytic Power')
%     hold on;
%     subplot(2,2,4)
    phase = angle(hilbert(GammaLFP(i,:)));
    fre = (phase(2:end)- phase(1:end-1))/1e-4/2/pi;
    plot(1e-4*t(2:end),fre(1e4+2:1e4+800))
    xlabel('Time(s)')
    ylabel('Analytic Frequency(rad/s)')
    hold on;
end
% for i = 1e4+1:tot
%     Power = flip(reshape(GLFPPower(:,i),[sqrt(no) sqrt(no)]));
%     surf(Power)
%     ts = sprintf('time = %8.1f ms', i*0.1);
%     title(ts);
%     pause(0.1)
% end
end