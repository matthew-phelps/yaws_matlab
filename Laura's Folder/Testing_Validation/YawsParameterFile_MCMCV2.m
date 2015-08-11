%% SEIR Model Parameters

% after the global statement, list all of your parameters.
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
        
        
       

rho = 0.043; %/12;
betaONE = .15;
betaTWO = betaONE;
mu = 0.0429;%;.016 remember to re-visit demographic params
alpha1 = .1295455;% see calculations in notebook pg1;
alpha2 = 1.234091;% see calculatins in notebook pg1;
gamma1 = 1/2.75;%see notebook   pg1         %1/9 - 1/2  % 3-6 weeks or more;
gamma2 = 1/6; % see notebook pg1
proportion = .095;
epsilon = 1.06510309;   %0.03617309
theta = 0.37598;  % see notebook pg2;           %1/6 - 1/0.5
psi =  (3.5/((33 - (0.905/epsilon))*35));
lambda =  9 * psi;         %1/20 - 1/2
eta =  2.5 * (psi + lambda);


% treatment in effect
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
phi1a = .9999999;
phi2a = .826;
phi3a =.826;
phi4a = .167;
phi1b = 1-phi1a;
phi2b = 1-phi2a;
phi3b = 1-phi3a;
phi4b = 1-phi4a;





    