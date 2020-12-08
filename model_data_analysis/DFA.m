function F = DFA(DATA,window)
% adjust from an online version
% detrended fluctuation analysis
N = length(DATA);
n = floor(N/window);
N1 = n*window;
% % overlap windows
% n = floor((N-window)/(0.5*window));
% N1 = (n+2)*0.5*window;

y = zeros(1,N1);
Yn = zeros(1,N1);
mean1 = mean(DATA(1:N1));
for i = 1:N1    
    y(i) = sum(DATA(1:i) - mean1);
end
for j = 1:n
    fitcoef = polyfit(1:window,y(((j-1)*window+1):j*window),1);    
    Yn(((j-1)*window+1):j*window) = polyval(fitcoef,1:window);
%     % overlap windows
%     fitcoef = polyfit(1:window,y(((j-1)*0.5*window+1):(j+1)*0.5*window),1);    
%     Yn(((j-1)*0.5*window+1):(j+1)*0.5*window) = polyval(fitcoef,1:window);
end
F = sqrt(mean(sum((y - Yn).^2)));
end