%% comparing the results of our model with the empirical results from
% Mitja in NEJM

% saving files to be read into R to make nice tables

% We need to show what percentage of people are still in the disease system
% after one round TCT at 82% coverage followed by one TTT at ~81% coverage


clear all
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results2x_TCTs_empircal_compare.mat')
results_mean = zeros(3,1);
std_result = zeros(3,1);
percent_pre_mean = zeros(3,1);
percent_pre_std = zeros(3,1);
for j = 1:3;
    %for l = 1:5;
    results_mean(j) = mean(results(:,j));
    std_result(j) = std(results(:,j));
    percent_pre_mean(j) = mean(percent_infect_pre(:,j));
    percent_pre_std(j) = std(percent_infect_pre(:,j));
    
    %end;
end;
hist(results(:,2))
save('empirical_compare_results_mean.mat', 'results_mean')
save('empirical_compare_std_result.mat', 'std_result')
save('empirical_compare_percent_pre_mean.mat', 'percent_pre_mean')
save('empirical_compare_percent_pre_std.mat', 'percent_pre_std')