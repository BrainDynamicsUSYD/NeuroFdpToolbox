% scatter plot based on reading from the txt.
fid = fopen( 'Scatter_Phi_EI_1.2_0.1_1.6.txt', 'r' );
cac = textscan( fid, '%s%n%s%n', 'Delimiter', ':' );
loop = cac{2};
scatter_val = cac{4};
loop_m = (find(sum(vec2mat(ismember([1:250],loop),10),2) == 10))';
scatter_vec = zeros(1,25);
for i = loop_m
    j = find(loop == (10*i - 9));
    scatter_vec(i) = mean(scatter_val(j:(j + 9)));
end
scatter_mat = vec2mat(scatter_vec,5);
figure(1)
imagesc(scatter_mat);
set(gca,'XTick',[1:5]);
set(gca,'XTickLabel',[1.2:0.1:1.6]);
xlabel('Phi_I');
set(gca,'YTick',[1:5]);
set(gca,'YTickLabel',[1.2:0.1:1.6]);
ylabel('Phi_E');
colorbar;
title('bump scatter');
saveas(gca,'bump_scatter.pdf');

% scatter_stdon = vec2mat(scatter_vec(1:2:49),5);
% scatter_stdoff = vec2mat(scatter_vec(2:2:50),5);
% figure(1)
% imagesc(scatter_stdon);
% set(gca,'XTick',[1.1:0.025:1.2]);
% set(gca,'YTick',[1.1:0.025:1.2]);
% colorbar;
% title('bump scatter STD on');
% saveas(gca,'bump_scatter_STD_on.pdf');
% figure(2)
% imagesc(scatter_stdoff);
% set(gca,'XTick',[1.1:0.025:1.2]);
% set(gca,'YTick',[1.1:0.025:1.2]);
% colorbar;
% title('bump scatter STD off');
% saveas(gca,'bump_scatter_STD_off.pdf');