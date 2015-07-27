%% Model File

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
        
%YawsParameterFile_TreatV2 %comment out when doing MCMC

startMonth = 1; 

timeBack = 400;

timeHorizon = 5000;%5000; 

timeSteps1 = zeros(1, timeBack+1);

timeSteps2 = zeros(1, timeHorizon+1);


for p=1:timeBack+1
    timeSteps1(p) = p-1;
end

for p=1:timeHorizon+1
    timeSteps2(p) = p-1;
end

%% initialize compartment values (using the names you assigned to your compartments in the ODE file

S_0 = 299; 

E_0 = 1;

IP1_0 = 0;

IP2_0 = 0;

L1_0 = 0;

IS_0 = 0;

L2_0 = 0;

IT_0 = 0;

R_0 = 0;




X_0 = zeros(1,8); 

X_0(1) = S_0; 

X_0(2) = E_0;

X_0(3) = IP1_0;

x_0(4) = IP2_0;

x_0(5) = L1_0;

X_0(6) = IS_0;

X_0(7) = L2_0;

X_0(8) = IT_0;

X_0(9) = R_0;

%% run ODE solver for timeSteps and initial values X_0

%[t2,X2]= ode15s('Yaws_ODE_FileV2', timeSteps2,X_0);
X21= ode15s(@Yaws_ODE_FileV2, timeSteps2,X_0);
X2 = deval(X21,timeSteps2)';
%% assign compartment names

X2 = max(0,X2);
%% assign compartment names (same as in ODE file)

S = X2(:,1); 
E = X2(:,2);
IP1 = X2(:,3);
IP2 = X2(:,4);
L1 = X2(:,5);
IS = X2(:,6);
L2 = X2(:,7);
IT = X2(:,8);
R = X2(:,9);

N = sum(X2(:,1:9),2); 
Nactive = X2(:,3) + X2(:,4) + X2(:,6);

%% version 1
% 
% ActivePrevPre = (X2(:,3)+X2(:,5))./N;
% ActivePrev = ActivePrevPre(end);

%% version 2: Oct 15

ActivePrevPre = (X2(:,3) + X2(:,4) + X2(:,6))./N;
ActivePrev = ActivePrevPre(end);

LatentPrevPre = (X2(:,5)+X2(:,7))./N;
LatentPrev = LatentPrevPre(end);

% LatentPrev/ActivePrev   %% comment out when doing MCMC
