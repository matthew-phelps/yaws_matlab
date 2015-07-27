%% FILES NOTES
% UNCERTAINTY analysis - looping model X number of times and sampling from posterior
% distribution of betaONE and Epsilon vect files AND the Efficacy of
% treatment all at the same time

% 1 treatment at all coverage levels and 744 months total time
% with "new" cure rate of 0.855

clear
tic
cbegin = fix(clock);
global  rho ...
    betaONE ...
    betaTWO ...
    mu ...
    alpha1 ...
    alpha2 ...
    gamma1 ...
    gamma2 ...
    proportion ...
    epsilon ...
    eta ...
    theta ...
    psi...
    lambda ...
    phi1a...
    phi2a...
    phi3a...
    phi4a...
    phi1b...
    phi2b...
    phi3b...
    phi4b...
    Ef...
    chi1...
    chi2...
    Rphi1a...
    Rphi2a...
    Rphi3a...
    Rphi4a...
    Rphi1b...
    Rphi2b...
    Rphi3b...
    Rphi4b...
    


YawsParameterFile_TreatV2
load('betaONEvect.mat')
load('epsilonVect.mat')


loops = 4200;
%%
% create vector of beta and epsilon values randomly drawn from MCMC vector
% output. We do this because we want the beta and epsilon to vary
% independtly during each run - so we can't use the original MCMC output, but we don't want
% to draw a random sample in each for-loop

% remove the burn-in period
betaONEvect_500 = betaONEvect(900:end);
epsilonvect_500 = epsilonvect(900:end);




%% generate vector of  of samples from parameters' FULL distributions
% epsilon & Beta found using LaplacesDemon package in R - p.interval(beta_500, HPD=F, MM=F, plot=T)
betaONEsample = datasample(betaONEvect_500, loops); % lower and upper bounds found using LaplacesDemon in R
epsilonSample = datasample(epsilonvect_500, loops); %epsilon cannot be below .0274 or you get spurious results
EfVect = binornd(250, 0.855, [loops,1])/250;
mean(EfVect);


coverage = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99];
coverage_names = {'Infect1x50', 'Infect1x55', 'Infect1x60', 'Infect1x65', 'Infect1x70','Infect1x75','Infect1x80','Infect1x85','Infect1x90','Infect1x95'};
results = zeros(loops, length(coverage));
people_infect = zeros(loops, length(coverage));
percent_infect_pre = zeros(loops, length(coverage));

%%

for j = 1:length(coverage) % iterate over the different coverage levels defined by the 'coverage variable
    
    for i = 1:loops
        betaONE = betaONEsample(i);  % sample randomly from Beta distribution
        epsilon = epsilonSample(i);  % sample randomly form Epsilon distribution
        psi =  (3.5/((33 - (0.905/epsilonSample(i)))*35));
        lambda =  9 * psi;         %1/20 - 1/2
        eta =  2.5 * (psi + lambda);
        
        
        startMonth = 1;
        
        timeHorizon = 744;
        
        timeSteps1 = [0:1:699];
        timeSteps2 = [0:1:1];
        timeSteps3 = [0:1:6];
        
        
        
        for k=1:length(timeSteps1)
            timeSteps1(k) = k-1;
        end
        
        for k=1:length(timeSteps2)
            timeSteps2(k) = k-1;
        end
        
        for k=1:length(timeSteps3)
            timeSteps3(k) = k-1;
        end
        
        
        %% initialize compartment values (using the names you assigned to your compartments in the ODE file
        
        S_0 =  267.2903;
        
        E_0 =  0.9955;
        
        IP1_0 = 1.0476;
        
        IP2_0 = 6.0476;
        
        L1_0 =  12.7306;
        
        IS_0 = 8.4223;
        
        L2_0 = 2.5995;
        
        IT_0 =  1.0075;
        
        R_0 =  0.9068;
        
        
        
        X_0 = zeros(1,9);
        
        X_0(1) = S_0;
        
        X_0(2) = E_0;
        
        X_0(3) = IP1_0;
        
        X_0(4) = IP2_0;
        
        X_0(5) = L1_0;
        
        X_0(6) = IS_0;
        
        X_0(7) = L2_0;
        
        X_0(8) = IT_0;
        
        X_0(9) = R_0;
        
        %% run ODE solver for timeSteps and initial values X_0
        
        
        % initial burn in
        timeSteps = [1:756];
        
        Ef = 0;
        Rphi1a = 0;
        Rphi2a = 0;
        Rphi3a = 0;
        Rphi4a = 0;
        Rphi1b = 0;
        Rphi2b = 0;
        Rphi3b = 0;
        Rphi4b = 0;
        chi1 = 0;
        chi2 = 0;
        [t1,X1]= ode15s('Yaws_ODE_FileV2',timeSteps1,X_0);
        
        %% 1st treatment
        
        Ef = 1;
        Rphi1a = -log(1-phi1a*EfVect(i)*coverage(j));
        Rphi2a = -log(1-phi2a*EfVect(i)*coverage(j));
        Rphi3a = -log(1-phi3a*EfVect(i)*coverage(j));
        Rphi4a = -log(1-phi4a*EfVect(i)*coverage(j));
        Rphi1b = -log(1-phi1b*EfVect(i)*coverage(j));
        Rphi2b = -log(1-phi2b*EfVect(i)*coverage(j));
        Rphi3b = -log(1-phi3b*EfVect(i)*coverage(j));
        Rphi4b = -log(1-phi4b*EfVect(i)*coverage(j));
        chi1 = 1;
        chi2 = 1;
        
        X21 = ode15s(@Yaws_ODE_FileV2,timeSteps2,X1(end,:));
        X2 = deval(X21,timeSteps2)';
        
        
        %% 3 year follow-up
        
        Ef = (0);
        Rphi1a = -log(1-phi1a);
        Rphi2a = -log(1-phi2a);
        Rphi3a = -log(1-phi3a);
        Rphi4a = -log(1-phi4a);
        Rphi1b = -log(1-phi1b);
        Rphi2b = -log(1-phi2b);
        Rphi3b = -log(1-phi3b);
        Rphi4b = -log(1-phi4b);
        chi1 = (0);
        chi2 = (0);
        
        X31 = ode15s(@Yaws_ODE_FileV2,timeSteps3,X2(end,:));
        X3 = deval(X31,timeSteps3)';
        
        
        X = [X1(1:end,:); X2(2:end,:); X3(2:end,:)];
        
        %% assign compartment names (same as in ODE file)
        
        S = X(:,1);
        E = X(:,2);
        IP1 = X(:,3);
        IP2 = X(:,4);
        L1 = X(:,5);
        IS = X(:,6);
        L2 = X(:,7);
        IT = X(:,8);
        R = X(:,9);
        
        N = sum(X(:,1:9),2);
        Nactive = X(:,3) + X(:,4) + X(:,6);
        Nlatent = X(:,5) + X(:,7);
        
        %% version 1 not tracking latent prev
        %
        % ActivePrevPre = (X2(:,3)+X2(:,5))./N;
        % ActivePrev = ActivePrevPre(end);
        
        %% version 2: Oct 15 - added tracking of latent prev
        
        ActivePrevPre = (X(:,3)+X(:,4)+X(:,6))./N;
        ActivePrev = ActivePrevPre(end);
        
        LatentPrevPre = (X(:,7)+X(:,5))./N;
        LatentPrev = LatentPrevPre(end);
        
        
        
        %% Sum number of people in system
        
        if E(end) < 0
            E(end) = 0;
        end;
        
        if IP1(end) < 0
            IP1(end) = 0;
        end;
        
        if IP2(end) < 0
            IP2(end) = 0;
        end;
        
        if L1(end) < 0
            L1(end) = 0;
        end;
        
        if L2(end) < 0
            L2(end) = 0;
        end;
        
        if IS(end) < 0
            IS(end) = 0;
        end;
        
        results(i,j) = ((E(end)+IP1(end)+ IP2(end) + L1(end) + L2(end) + IS(end)) / N(end));
        people_infect(i,j) = ((E(end)+IP1(end)+ IP2(end) + L1(end) + L2(end) + IS(end)) );
        percent_infect_pre (i, j) = ((E(699) + IP1(699) + IP2(699) + L1(699) + L2(699) + IS(699)) / N(699));
    end;
    
end;




% %% Prevalence reduction
% 
% prev_post_mean = zeros(1, length(coverage));
% prev_pre_mean = zeros(1, length(coverage));
% ts_lower_post = zeros(1, length(coverage));
% 
% CI_lower_post = zeros(1, length(coverage));
% ts_upper_post = zeros(1, length(coverage));
% CI_upper_post = zeros(1, length(coverage));
% SEM_pre = zeros(1, length(coverage));
% ts_lower_pre = zeros(1, length(coverage));
% ts_upper_pre = zeros(1, length(coverage));
% CI_lower_pre = zeros(1, length(coverage));
% CI_upper_pre = zeros(1, length(coverage));
% SEM_post = zeros(1, length(coverage));
% mean_diff = zeros(1, length(coverage));
% max_diff = zeros(1, length(coverage));
% min_diff = zeros(1, length(coverage));
% proportion_reduction_pe = zeros(1, length(coverage));
% proportion_reduction_MAX = zeros(1, length(coverage));
% proportion_reduction_MIN = zeros(1, length(coverage));
% 
% for j = 1:length(coverage);
%     
%     prev_pre_mean(j) = mean(percent_infect_pre(:,j));
%     prev_post_mean(j) = mean(results(:,j));
%     
%     SEM_pre(j) = std(percent_infect_pre(:,j)) / sqrt(loops);
%     ts_lower_pre(j) = tinv([0.025], loops-1);
%     CI_lower_pre(j) = prev_pre_mean(j) + ts_lower_pre(j) * SEM_pre(j);
%     ts_upper_pre(j) = tinv([0.975], loops-1);
%     CI_upper_pre(j) = prev_pre_mean(j) + ts_upper_pre(j) * SEM_pre(j);
%     
%     SEM_post(j) = std(results(:,j)) / sqrt(loops);
%     ts_lower_post(j) = tinv([0.025], loops-1);
%     CI_lower_post(j) = prev_post_mean(j) + ts_lower_post(j) * SEM_post(j);
%     ts_upper_post(j) = tinv([0.975], loops-1);
%     CI_upper_post(j) = prev_post_mean(j) + ts_upper_post(j) * SEM_post(j);
%     
%     mean_diff(j) = prev_pre_mean(j) - prev_post_mean(j);
%     max_diff(j) = CI_upper_pre(j) - CI_lower_post(j);
%     min_diff(j) = CI_lower_pre(j) - CI_upper_post(j);
%     
%     proportion_reduction_pe(j) = mean_diff(j) / prev_pre_mean(j);
%     proportion_reduction_MAX(j) = max_diff(j) / prev_pre_mean(j);
%     proportion_reduction_MIN(j) = min_diff(j) / prev_pre_mean(j);
% end;
% 
% 
% %% Saving output
% save('x1_TCT_prev_reduc_mean.mat', 'proportion_reduction_pe')
% save('x1_TCT_prev_reduc_MAX.mat', 'proportion_reduction_MAX')
% save('x1_TCT_prev_reduc_MIN.mat', 'proportion_reduction_MIN')

%% prob Success
success = zeros(loops, length(coverage));
for j = 1:length(coverage);
    for i = 1:loops;
        if people_infect(i,j) < 0.5;
            success(i,j) = 1;
        else success(i,j) = 0;
        end;
    end;
end;

%% prob success

prob_success = zeros(1, length(coverage));
for j = 1:length(coverage);
    prob_success(1,j) = sum(success(:,j)) / loops;
end;

save results1x_TCTs_empirical_compare.mat

toc