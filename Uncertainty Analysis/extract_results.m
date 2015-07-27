clear all
%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results3x_TCTs.mat');
save('x3_TCT_Results.mat', 'prob_success')
clear all
%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results4x_TCTs.mat');
save('x4_TCT_Results.mat', 'prob_success')

clear all
%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results5x_TCTs_1p.mat');
save('x5_TCT_Results_1p.mat', 'prob_success')

clear all
%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results6x_TCTs.mat');
save('x6_TCT_Results.mat', 'prob_success')



%%
%%





%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results3x_Pr.mat');
people_infect_avg = zeros(5,10);
for i = 1:10;
    for j = 1:5;
        people_infect_avg(j,i) = mean(people_infect(:,i,j));
    end;
end;
save('x3_TTT_ppl.mat', 'people_infect_avg')
clear all

%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results4x_varTime.mat');
people_infect_avg = zeros(5,10);
for i = 1:10;
    for j = 1:5;
        people_infect_avg(j,i) = mean(people_infect(:,i,j));
    end;
end;
save('x4_TTT_ppl.mat', 'people_infect_avg')
clear all

%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results5x_varTime.mat');
people_infect_avg = zeros(5,10);
for i = 1:10;
    for j = 1:5;
        people_infect_avg(j,i) = mean(people_infect(:,i,j));
    end;
end;
save('x5_TTT_ppl.mat', 'people_infect_avg')
clear all


%%
%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results6x_varTime.mat');
people_infect_avg = zeros(5,10);
for i = 1:10;
    for j = 1:5;
        people_infect_avg(j,i) = mean(people_infect(:,i,j));
    end;
end;
save('x6_TTT_ppl.mat', 'people_infect_avg')
clear all






%%
load('C:\Users\wrz741\Google Drive\Yaws project\MATLAB\Uncertainty Analysis\results6x_hybrid.mat');
save('x6_hybrid_Results.mat', 'prob_success')