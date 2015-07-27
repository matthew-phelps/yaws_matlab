%% SEIR Model ODE Equations

function Xdot = Yaws_ODE_File(t,X) %use your file name

%
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

%% assign compartment names

% These refer to my compartments (Scu80 - susceptible, unvaccinated in the
% low risk category); feel free to be as descriptive as you want with yours

S = X(1);
E = X(2);
IP = X(3);
L1 = X(4);
IS = X(5);
L2 = X(6);
IT = X(7);
R = X(8);

%% ODEs

Xdot = zeros(size(X));  
Ntot = S + E + IP + L1 + IS + L2 + IT + R;

Force = (betaONE*IP+betaTWO*IS)/Ntot;

% Differential equations:
%dS/dt
Xdot(1) = (rho*Ntot) - (S*Force) - (mu*S);
 
%dE/dt 
Xdot(2) = (S*Force) - (alpha*E) - (mu*E);
 
%dIP/dt 
Xdot(3) = (alpha*E) - (gamma*IP) - (mu*IP);

%dL1/dt
Xdot(4) = ((1-proportion)*gamma*IP) - epsilon*L1 - mu*L1;

%dIS/dt
Xdot(5) = eta*L2 -theta*IS + proportion*gamma*IP + epsilon*L1 - mu*IS;
 
%dL2/dt
Xdot(6) = theta*IS - 2*lambda*L2 - eta*L2 - mu*L2;
 
%dIT/dt
Xdot(7) = lambda*L2 - mu*IT;

%dR/dt
Xdot(8) = lambda*L2 - mu*R;
 
end