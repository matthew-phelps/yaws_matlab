%% SEIR Model Parameters

% after the global statement, list all of your parameters.
global  rho ...
        betaONE ...
        betaTWO ... 
        mu ...
        alpha ...
        gamma ...
        proportion ...
        epsilon ...
        eta ...
        theta ...
        lambda ...       

rho = 0.043; %/12;
betaONE = 0.8;
betaTWO = betaONE;
mu = 0.04;%;.016 remember to re-visit demographic params
alpha = .7;%/12;
gamma = 1/.75;            %1/9 - 1/3  % 3-6 weeks or more;
proportion = .095;
epsilon = 1/5; %1/100;         %1/24 - 0
eta = 0.1;             %1/60-1/3
theta = 1/3;           %1/6 - 1/0.5
lambda =  1/5;           %1/20 - 1/2


    