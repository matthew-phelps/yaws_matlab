%% Model File

% Copy same global list that you included in your ODE and Parameter files
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
        
YawsParameterFile_aug % This will need to change if you rename the parameter file

% Just keep the next 11 lines of code and I will explain what they do when
% we meet next
startMonth = 1; 

timeBack = 400;

timeHorizon = 5000; 

timeSteps1 = zeros(1, timeBack+1);

timeSteps2 = zeros(1, timeHorizon+1);

ActivePrev = zeros(timeHorizon,1); %creating empty matrix to be populated by prevelance values for the given timesteps

for p=1:timeBack+1
    timeSteps1(p) = p-1;
end

for p=1:timeHorizon+1
    timeSteps2(p) = p-1;
end

%% initialize compartment values (using the names you assigned to your compartments in the ODE file

S_0 = 300; 

E_0 = 1;

IP_0 = 0;

L1_0 = 0;

IS_0 = 1;

L2_0 = 0;

IT_0 = 0;

R_0 = 0;



X_0 = zeros(1,8); 

X_0(1) = S_0; 

X_0(2) = E_0;

X_0(3) = IP_0;

x_0(4) = L1_0;

X_0(5) = IS_0;

X_0(6) = L2_0;

X_0(7) = IT_0;

X_0(8) = R_0;

%% run ODE solver for timeSteps and initial values X_0

[t2,X2]= ode45('Yaws_ODE_File', timeSteps2,X_0);

%% assign compartment names (same as in ODE file)

S = X2(:,1); % 1st column corresponds to Unvaccinated Susceptible Children
E = X2(:,2);
IP = X2(:,3);
L1 = X2(:,4);
IS = X2(:,5);
L2 = X2(:,6);
IT = X2(:,7);
R = X2(:,8);


N = X2(:,1) + X2(:,2) + X2(:,3) + X2(:,4) + X2(:,5) + X2(:,6)+ X2(:,7) + X2(:,8); % you only need to include this if you want tokeep track of total population size

% defining prevelance
for j=1:timeHorizon; %timeHorizon = the end, the point at which you reach an equilibrium, user defined value
   ActivePrev(j) = (IS(end)+ IP(end))/N(end);
end

%ChildrenPrev(end);

 
% If you want to plot individual compartments over time, you can use the
% following plot statement
 
% figure(1)
% plot(timeSteps2,IP,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('IP')
% 
% figure(2)
% plot(timeSteps2,IS,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('IS')
% 
% figure(3)
% plot(timeSteps2,L,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('L')
% 
% figure(4)
% plot(timeSteps2,R,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('R')
% 
% figure(5)
% plot(timeSteps2,N,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('N')
% 
% figure(6)
% plot(timeSteps2,S,'-b') %where Scu will be replaced with one of your compartment names
% xlabel('Month')
% ylabel('Number of Individuals')
% legend('S')