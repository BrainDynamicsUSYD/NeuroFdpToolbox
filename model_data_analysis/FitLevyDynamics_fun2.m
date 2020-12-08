
function [distance,delta_x,factor] = FitLevyDynamics_fun2(R,mode,bin,jump_win)
close all;
% clc;
% mode = 'bayesian';


hw = 31;
fw = 2*hw+1;

% donot use bayes since it cannot generate accurate jump size in a fixed
% window
switch mode
    case 'bayes'
        %         if 1
        if ~isfield(R,'grid') || ~isfield(R.grid,'bayes') ||...
                ~isfield(R.grid.bayes,'jump_win') ||...
                R.grid.bayes.jump_win ~= jump_win ||...
                R.grid.bayes.bin_size ~= bin
            R = get_grid_firing_centreYL(R,'mode','bayesian','win_len',bin,'jump_win',jump_win/10);
            SaveRYG_nocell(R)
        end
%         % no gap jumping
%         % get the jump size and distance: when calculating these, position vector
%         % must belong to one complete trajectory, i.e. no NaN in them
%         [~, ~, du,~,start] = seq_postprocess(isnan(R.grid.quick.radius),1);
%         pos = R.grid.quick.centre'/hw*pi;
%         jump_win = jump_win/10; % ms
%         
%         jump_size = [];
%         for i = 1:length(du)
%             p = pos(start(i):(start(i)+du(i)-1),:);
%             jump_s = wrapToPi(p(1:end-jump_win,:)-p((jump_win+1):end,:));
% %             if strcmp(mode,'bayes')
% %                 factor = R.grid.bayes.bayes_factor_ln(start(i):(start(i)+du(i)-1));
% %                 factor(factor <= log(100)) = 0;
% %                 factor(factor > 0) = 1;
% %                 ff = factor(1:end-jump_win).*factor((jump_win+1):end);
% %                 jump_s = jump_s(ff>0,:);
% %             end
%             jump_size = [jump_size;jump_s];
%         end
%         
%         jump_size = wrapToPi(pos(1:end-jump_win,:)-pos((jump_win+1):end,:));
%         jump_dist = sqrt(sum(jump_size.*jump_size,2)); %radial distance of increment
%         % convert back to real coordinate
%         R.grid.bayes.jump_size = jump_size(:).*fw/2/pi;
%         R.grid.bayes.jump_dist = jump_dist(:).*fw/2/pi;
        
        distance = R.grid.bayes.jump_dist;
        delta_x =R.grid.bayes.jump_size;
        factor = R.grid.bayes.bayes_factor_ln;
    case 'quick'
%                 if 1
        if ~isfield(R,'grid') || ~isfield(R.grid,'quick') ||...
                ~isfield(R.grid.quick,'jump_win') ||...
                R.grid.quick.jump_win ~= jump_win ||...
                R.grid.quick.bin_size ~= bin
            R = get_grid_firing_centreYL(R,'mode','quick','win_len',bin,'jump_win',jump_win/10);
            SaveRYG_nocell(R)
        end
%         % no gap jumping
%         % get the jump size and distance: when calculating these, position vector
%         % must belong to one complete trajectory, i.e. no NaN in them
%         [~, ~, du,~,start] = seq_postprocess(isnan(R.grid.quick.radius),1);
%         pos = R.grid.quick.centre'/hw*pi;
%         jump_win = jump_win/10; % ms        
%         
%         jump_size = [];
%         for i = 1:length(du)
%             p = pos(start(i):(start(i)+du(i)-1),:);
%             jump_s = wrapToPi(p(1:end-jump_win,:)-p((jump_win+1):end,:));
%             jump_size = [jump_size;jump_s];
%         end
%         
%         jump_size = wrapToPi(pos(1:end-jump_win,:)-pos((jump_win+1):end,:));
%         jump_dist = sqrt(sum(jump_size.*jump_size,2)); %radial distance of increment
%         % convert back to real coordinate
%         R.grid.quick.jump_size = jump_size(:).*fw/2/pi;
%         R.grid.quick.jump_dist = jump_dist(:).*fw/2/pi;
%         
        distance = R.grid.quick.jump_dist;
        delta_x =R.grid.quick.jump_size;
        factor = [];
    case 'pop'
        if ~isfield(R.grid,'population_vector_centre') ||...
                ~isfield(R.grid.population_vector_centre,'jump_win') ||...
                R.grid.population_vector_centre.jump_win ~= jump_win/10 ||...
                R.grid.population_vector_centre.bin_size ~= bin/10
            R = get_grid_firing_centre_populationvector(R,'tau',bin/10,'jump_win',jump_win/10);
            %             SaveRYG_nocell(R)
        end
        distance = R.grid.population_vector_centre.jump_dist;
        delta_x = R.grid.population_vector_centre.jump_size;
        factor = [];
    case 'real-time_I'
        x = R.neuron_stats.pos_x_I{1};
        y = R.neuron_stats.pos_y_I{1};
        
        distance = wrapToPi((diff(x).^2+diff(y).^2).^0.5).*63/2/pi;
        delta_x = [wrapToPi(diff(x));wrapToPi(diff(y))].*63/2/pi;
    case 'real-time_V'
        x = R.neuron_stats.pos_x_V{1};
        y = R.neuron_stats.pos_y_V{1};
        
        distance = wrapToPi((diff(x).^2+diff(y).^2).^0.5).*63/2/pi;
        delta_x = [wrapToPi(diff(x));wrapToPi(diff(y))].*63/2/pi;
    otherwise
        error('Please input correct mode!')
end


% delta_x = delta_x(:);
% if increment_flag %if true, use fit increment to SaS distribution
%     pd_incre = sasML(delta_x,'sas');
%     pd_dist = sasML(distance,'asas');
% else
%     pd_incre = sasML(delta_x,'stable');
%     pd_dist = sasML(distance,'stable');
% end
% x_variable = unique(delta_x);
% y_variable = pdf(pd_incre,x_variable);
% kurto_dis = kurtosis(distance,0);
% skw_dis = skewness(distance,0);
% mean_dis = mean(distance,'omitnan');
% figure('visible','on')
% plot(x_variable,y_variable)
% hold on
%
% % plot(edge(1:end-1),N)
% edge_custom =linspace(-pi,pi,31);
% h=histogram(delta_x,edge_custom,'normalization','pdf');% #bins:60
% xlabel('X cahnge (grid points)')
% ylabel('Probability')
% legend({'fitted','data'})
% title('X change')
% saveas(gcf,[num2str(bin),'_',datestr(now,'yyyymmddHHMM-SSFFF'),'.jpg'])
%
% figure('visible','on')
% edge_custom =linspace(0,max(distance),30);
% [N,edge]=histcounts(distance,edge_custom,'normalization','pdf');
% loglog(edge(1:end-1),N)
% title('Distance')
% xlabel('Distance (grid points)')
% ylabel('Probability')
% % histogram(distance,'normalization','probability');
% saveas(gcf,[num2str(bin),'_',datestr(now,'yyyymmddHHMM-SSFFF'),'.jpg'])
end

function SaveRYG_nocell( Result )
%This function saves Result_cell as explanatory and response
%variables to .mat file
%   Detailed explanation goes here

name_temp = strcat( Result.stamp, '_RYG.mat');
fprintf('Saving results into file %s...\n', name_temp);
save(name_temp, '-struct', 'Result', '-v7.3'); % -v7.3 for >2GB

disp('Saving done.')

end