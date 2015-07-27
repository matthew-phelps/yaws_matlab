%% Get 95% distribution of the two fitted parameters and Treatment Efficacy
clear all
close all

load('betaONEvect.mat')
load('epsilonVect.mat')

betaONEvect_500 = betaONEvect(500:end);
epsilonvect_500 = epsilonvect(500:end);

%histfit(betaONEvect_500, 500, 'normal')
%histfit(epsilonvect_500, 500, 'lognormal')
% plot (epsilonvect(500:end))
% plot(betaONEvect(150:end))

% fit distributions to data
beta_normal = fitdist(betaONEvect_500, 'Normal');
epsilon_lognrm = fitdist(epsilonvect_500, 'lognormal');
epsilon_kernel = fitdist(epsilonvect_500, 'kernel');

% plot distributuions
% 
 x = min(betaONEvect_500):0.001:max(betaONEvect_500);
 y = pdf(beta_normal, x);
% 
plot(x,y)

% generate vector of 100k of samples from beta's distribution
beta_sample_norm = random(beta_normal, 100000,1);


% generate vector of 100k of samples from epsilon's distribution
epsilon_sample_lognrm = random(epsilon_lognrm, 100000,1);
epsilon_sample_kernel = random(epsilon_kernel, 100000,1);


% save vectors to file for input into UA

save sample_vectors epsilon_sample_lognrm epsilon_sample_kernel beta_sample_norm

[muhat, sigmahat] = normfit (betaONEvect_500)
muhat_hi = muhat + 1.96*sigmahat
muhat_low = muhat - 1.96*sigmahat

posterior_moments(betaONEvect_500, 1, 0.95)





