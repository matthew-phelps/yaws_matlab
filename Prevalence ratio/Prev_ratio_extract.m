%% Prevalence reduction
clear
load('results1x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(106:4095,:);
prev_post_95 = sorted_prev_post(106:4095,:);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = min(prev_ratio_vect);
prev_ratio_max = max(prev_ratio_vect);
prev_ratio_PE = mean(prev_ratio_vect);
%%


%%
save('x1_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x1_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x1_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')






%% Prevalence reduction
clear
load('results2x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
% %%
%%
save('x2_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x2_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x2_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results3x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
%%
save('x3_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x3_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x3_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results4x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%

%%
save('x4_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x4_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x4_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results5x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x5_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x5_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x5_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results6x_TCTs_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x6_TCT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x6_TCT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x6_TCT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results2x_TTT_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x2_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x2_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x2_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results3x_TTT_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x3_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x3_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x3_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results4x_TTT_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x4_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x4_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x4_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results5x_TTT_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x5_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x5_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x5_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')



%%
%%
%% Prevalence reduction
clear
load('results6x_TTT_empirical_compare.mat')
sorted_prev_pre = sort(percent_infect_pre);
sorted_prev_post = sort(results);
prev_pre_95 = sorted_prev_pre(31:1170, :, :);
prev_post_95 = sorted_prev_post(31:1170, :, :);

prev_ratio_vect = prev_post_95 ./ prev_pre_95;
prev_ratio_min = (reshape(min(prev_ratio_vect), [10,5])).';
prev_ratio_max = (reshape(max(prev_ratio_vect), [10,5])).';
prev_ratio_PE = (reshape(mean(prev_ratio_vect), [10,5])).';
%%
save('x6_TTT_prev_ratio_PE.mat', 'prev_ratio_PE')
save('x6_TTT_prev_ratio_MAX.mat', 'prev_ratio_max')
save('x6_TTT_prev_ratio_MIN.mat', 'prev_ratio_min')