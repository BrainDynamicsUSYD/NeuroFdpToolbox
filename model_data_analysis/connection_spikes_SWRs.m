% check the connection between SWRs with spikes' number in the same time
[is_SWR,~] = CollectVectorYG('LFP','LFP.ripple_event.is_SWR');
is = vec2mat(is_SWR,16);
[spikes_E,~] = CollectVectorYG('num_spikes','num_spikes{1}');
[spikes_I,~] = CollectVectorYG('num_spikes','num_spikes{2}');
sum_SWR = sum(is,2);
% is_SWR = R.LFP.ripple_event.is_SWR;
% spikes_E = R.num_spikes{1};
% spikes_I = R.num_spikes{2};
% sum_SWR = sum(is_SWR);
I = [];
J = [];
ave_SWR = [];
n_spikes_E = [];
n_spikes_I = [];
for t = 1:(length(sum_SWR) -1)
    if sum_SWR(t) == 0 && sum_SWR(t + 1) > 0
        I = [I t + 1];
    elseif sum_SWR(t) > 0 && sum_SWR(t + 1) == 0
        J = [J t];
    end
end
for t = 1:length(I)
    n_spikes_E = [n_spikes_E sum(spikes_E(I(t):J(t)))];
    n_spikes_I = [n_spikes_I sum(spikes_I(I(t):J(t)))];
    ave_SWR = [ave_SWR mean(sum_SWR(I(t):J(t)))];
end
n_spikes_EI = n_spikes_E + n_spikes_I;
nn_spikes_E = round(n_spikes_E./ave_SWR);
nn_spikes_I = round(n_spikes_I./ave_SWR);
nn_spikes_EI = nn_spikes_E + nn_spikes_I;
figure(1)
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Connection between Spikes and SWRs', 'FontSize', 12, 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
subplot(2,3,1)
histogram(n_spikes_E,25)
xlabel({'Excitatory',['min=',num2str(min(n_spikes_E)),',mean=',num2str(round(mean(n_spikes_E)))]},'FontSize',10)
subplot(2,3,2)
histogram(n_spikes_I,25)
xlabel({'EInhibitory',['min=',num2str(min(n_spikes_I)),',mean=',num2str(round(mean(n_spikes_I)))]},'FontSize',10)
subplot(2,3,3)
histogram(n_spikes_EI,25)
xlabel({'All',['min=',num2str(min(n_spikes_EI)),',mean=',num2str(round(mean(n_spikes_EI)))]},'FontSize',10)
subplot(2,3,4)
histogram(nn_spikes_E,25)
xlabel({'Excitatory(normalized)',['min=',num2str(min(nn_spikes_E)),',mean=',num2str(round(mean(nn_spikes_E)))]},'FontSize',10)
subplot(2,3,5)
histogram(nn_spikes_I,25)
xlabel({'Inhibitory(normalized)',['min=',num2str(min(nn_spikes_I)),',mean=',num2str(round(mean(nn_spikes_I)))]},'FontSize',10)
subplot(2,3,6)
histogram(nn_spikes_EI,25)
xlabel({'All(normalized)',['min=',num2str(min(nn_spikes_EI)),',mean=',num2str(round(mean(nn_spikes_EI)))]},'FontSize',10)