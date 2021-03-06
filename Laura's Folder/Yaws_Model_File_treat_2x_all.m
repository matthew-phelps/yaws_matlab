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



coverage = [0.50,0.55,0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95];
coverage_names = {'Infect1x50', 'Infect1x55', 'Infect1x60', 'Infect1x65', 'Infect1x70','Infect1x75','Infect1x80','Infect1x85','Infect1x90','Infect1x95'}; 
PE_2x = zeros(length(coverage),1);



for j = 1:length(coverage) % iterate over the different coverage levels defined by the 'coverage variable

startMonth = 1; 



timeHorizon = 744; 


% 
% timeSteps1 = zeros(1,700);
% 
% timeSteps2 = zeros(1,2);
% 
% timeSteps3 = zeros(1,7);
% 
% timeSteps4 = zeros(1,2);
% 
% timeSteps5 = zeros(1,37);

timeSteps1 = [0:1:699];
timeSteps2 = [0:1:1];
timeSteps3 = [0:1:6];
timeSteps4 = [0:1:1];
timeSteps5 = [0:1:36];



for k=1:length(timeSteps1)
    timeSteps1(k) = k-1;
end

for k=1:length(timeSteps2)
    timeSteps2(k) = k-1;
end

for k=1:length(timeSteps3)
    timeSteps3(k) = k-1;
end

for k=1:length(timeSteps4)
    timeSteps4(k) = k-1;
end

for k=1:length(timeSteps5)
    timeSteps5(k) = k-1;
end

%end
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



X_0 = zeros(1,8); 

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
timeSteps = [1:744]; 

% initial burn in

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

Ef = -log(1-0.855);
Rphi1a = -log(1-phi1a);
Rphi2a = -log(1-phi2a);
Rphi3a = -log(1-phi3a);
Rphi4a = -log(1-phi4a);
Rphi1b = -log(1-phi1b);
Rphi2b = -log(1-phi2b);
Rphi3b = -log(1-phi3b);
Rphi4b = -log(1-phi4b);
chi1 = -log(1-coverage(j));
chi2 = -log(1-coverage(j));

X21 = ode15s(@Yaws_ODE_FileV2,timeSteps2,X1(end,:));
X2 = deval(X21,timeSteps2)';

%% 6month inter-period

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

%% Treatment #2
% 
Ef = -log(1-.855);
Rphi1a = -log(1-phi1a);
Rphi2a = -log(1-phi2a);
Rphi3a = -log(1-phi3a);
Rphi4a = -log(1-phi4a);
Rphi1b = -log(1-phi1b);
Rphi2b = -log(1-phi2b);
Rphi3b = -log(1-phi3b);
Rphi4b = -log(1-phi4b);
chi1 = -log(1-coverage(j));
chi2 = (0);

X41 = ode15s(@Yaws_ODE_FileV2,timeSteps4,X3(end,:));
X4 = deval(X41, timeSteps4)';
%% 3 year follow-up

Ef = 0;
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

[t5,X5]= ode15s('Yaws_ODE_FileV2',timeSteps5,X4(end,:));

X = [X1(1:end,:); X2(2:end,:); X3(2:end,:); X4(2:end,:); X5(2:end,:)];






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


 PE_2x(j) = (E(end)+IP1(end)+ IP2(end) + L1(end) + L2(end) + IS(end))/ N(end);





end;
save pe_2x.mat PE_2x

toc