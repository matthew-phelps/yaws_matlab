clear
%% Prevalence reduction
load('results3x_TCTs_empirical_compare.mat')
%%
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :); % remove the first 30 and last 30 obs
prev_post_95 = sorted_prev_post(31:1170, :, :);

prop_reduction_vec = (percent_infect_pre - results)./percent_infect_pre;
proportion_reduction_MIN = min(prop_reduction_vec);
proportion_reduction_MAX = max(prop_reduction_vec);
proportion_reduction_pe = mean(prop_reduction_vec);

%%
prev_post_mean = zeros(length(time),length(coverage));
prev_pre_mean = zeros(length(time),length(coverage));
CI_lower_post = zeros(length(time),length(coverage));
CI_upper_post = zeros(length(time),length(coverage));
CI_lower_pre = zeros(length(time),length(coverage));
CI_upper_pre = zeros(length(time),length(coverage));
prev_ratio_PE = zeros(length(time),length(coverage));
prev_ratio_min = zeros(length(time),length(coverage));
prev_ratio_max = zeros(length(time),length(coverage));

for j = 1:length(coverage);
    for l = 1:length(time);
        
        prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
        prev_post_mean(l,j) = mean(results(:,j,l));
        
        CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
        CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
        
        CI_lower_post(l,j) = min(prev_post_95(:, j,l));
        CI_upper_post(l,j) = max(prev_post_95(:, j,l));
        
        prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
        prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
        prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
        
        mean_diff(l,j) = prev_pre_mean(l,j) - prev_post_mean(l,j);
        max_diff(l,j) = CI_upper_pre(l,j) - CI_lower_post(l,j);
        min_diff(l,j) = CI_lower_pre(l,j) - CI_upper_post(l,j);
        
        proportion_reduction_pe(l,j) = mean_diff(l,j) / prev_pre_mean(l,j);
        proportion_reduction_MAX(l,j) = max_diff(l,j) / prev_pre_mean(l,j);
        proportion_reduction_MIN(l,j) = min_diff(l,j) / prev_pre_mean(l,j);

    end;
end;

%%
% 
% 
% 
% 
% %%
% %%
% %% Prevalence reduction
% load('results3x_TCTs_empirical_compare.mat')
% 
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x3_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x3_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x3_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %%
% %% Prevalence reduction
% load('results4x_TCTs_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x4_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x4_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x4_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %%
% %% Prevalence reduction
% load('results5x_TCTs_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x5_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x5_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x5_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %%
% %% Prevalence reduction
% load('results6x_TCTs_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x6_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x6_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x6_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %%
% %% Prevalence reduction
% load('results2x_TTT_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x2_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x2_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x2_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %% 
% %% Prevalence reduction
% clear
% load('results3x_TTT_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x3_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x3_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x3_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %% 
% %% Prevalence reduction
% clear
% load('results4x_TTT_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x4_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x4_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x4_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %% 
% %% Prevalence reduction
% clear
% load('results5x_TTT_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x5_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x5_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x5_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')
% 
% 
% 
% %%
% %% 
% %% Prevalence reduction
% clear
% load('results6x_TTT_empirical_compare.mat')
% sorted_prev_pre = sort(percent_infect_pre);
% sorted_prev_post = sort(results);
% prev_pre_95 = sorted_prev_pre(31:1170, :, :);
% prev_post_95 = sorted_prev_post(31:1170, :, :);
% 
% 
% %%
% prev_post_mean = zeros(length(time),length(coverage));
% prev_pre_mean = zeros(length(time),length(coverage));
% CI_lower_post = zeros(length(time),length(coverage));
% CI_upper_post = zeros(length(time),length(coverage));
% CI_lower_pre = zeros(length(time),length(coverage));
% CI_upper_pre = zeros(length(time),length(coverage));
% prev_ratio_PE = zeros(length(time),length(coverage));
% prev_ratio_min = zeros(length(time),length(coverage));
% prev_ratio_max = zeros(length(time),length(coverage));
% 
% for j = 1:length(coverage);
%     for l = 1:length(time);
%         
%         prev_pre_mean(l,j) = mean(percent_infect_pre(:,j,l));
%         prev_post_mean(l,j) = mean(results(:,j,l));
%         
%         CI_lower_pre(l,j) = min(prev_pre_95(:,j,l));
%         CI_upper_pre(l,j) = max(prev_pre_95(:, j,l));
%         
%         CI_lower_post(l,j) = min(prev_post_95(:, j,l));
%         CI_upper_post(l,j) = max(prev_post_95(:, j,l));
%         
%         prev_ratio_min(l,j) = CI_lower_post(l,j) / CI_upper_pre(l,j);
%         prev_ratio_PE(l,j) = prev_post_mean(l,j) / prev_pre_mean(l,j);
%         prev_ratio_max(l,j) = CI_upper_post(l,j) / CI_lower_pre(l,j);
% 
%     end;
% end;
% 
% %%
% save('x6_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
% save('x6_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
% save('x6_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')