format compact
tic

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
        
loops2 = 100000 ; %200150;

ActivePrevvec = zeros(loops2+1,1); %creating space to be populated with the prevalance values from the model. Loops is in reference to how many iterations we will be doing

LatentPrevvec = zeros(loops2+1,1); % remove if we don't use ratio method 

LactivePrev = zeros(loops2+1,1); %This is space where likelihood values will be stored
LlatentPrev = zeros(loops2+1,1);

Ltot = zeros(loops2+1,1); % this is the likelihood of both active and latent combined

YawsParameterFile_MCMCV2; % my parameter file

Yaws_Model_File_baseV2; %my model file

ActivePrevvec(1) = ActivePrev(end); %we only save values that are accepted by the Hastings algorithm. 

LatentPrevvec(1) = LatentPrev(end); % remove if we don't use ratio method

%% These are the three parameters we want to fit

betaONEvect = zeros(loops2+1,1);
betaONEvect(1) = betaONE;      % 0.000001 - 0.00001
move_betaONE = zeros(loops2,1);
step_betaONE = 0.08*betaONEvect(1);

epsilonvect = zeros(loops2+1,1);
epsilonvect(1) = epsilon;      % 0.000001 - 0.00001
move_epsilon = zeros(loops2,1);
step_epsilon = 0.08*epsilonvect(1);

%% LIKELIHOOD FUNCTIONS

LactivePrev_old = log(binopdf(round(13490*ActivePrevvec(1)),13490,0.024));% the probability (P) comes from the data, the x comes from the model, N comes
% from the study population from which P was estimated
% N = 1176; P=.047
LlatentPrev_old = log(betapdf((ActivePrevvec(1)/LatentPrevvec(1)),9.78,19.56)); % To restrain so that 95% of the time the ratio remains between 2 to 6
Ltot_old = LactivePrev_old + LlatentPrev_old;


%Take the x output from the model and use that as the new input
%(round(1176*prev), 1776, 0.047)

LactivePrev(1) = LactivePrev_old;
LlatentPrev(1) = LlatentPrev_old;
Ltot(1) = Ltot_old;

Aa = [1:2:700150];%the big number at the end is the total number of loops we do
Bb = [2:2:700150];

%% For loops


for i = 1:loops2 
     %% BetaONE   
    
    if find(Aa ==i)>0

        triang1 = makedist('Triangular','a',betaONEvect(i)-step_betaONE,'b',betaONEvect(i),'c',betaONEvect(i)+step_betaONE);
        %rng('default'); % For reproducibility
        betaONEvect(i+1) = random(triang1,1,1);

        %betaONEvect(i+1) = betaONEvect(i) - step_betaONE + rand(1)*2*step_betaONE;      
         
        epsilonvect(i+1)    = epsilonvect(i);
        
        if betaONEvect(i+1) <= 0 || betaONEvect(i+1) > 50; %0.0002, 0.002

            betaONEvect(i+1) = betaONEvect(i);
            
            move_betaONE(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            LactivePrev(i+1) = LactivePrev_old;
            LlatentPrev(i+1) = LlatentPrev_old;
            Ltot(i+1) = Ltot_old;
            
        else 
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            epsilon = epsilonvect(i+1);
            psi =  (3.5/((33 - (0.905/epsilon))*35));
            lambda =  9 * psi;         %1/20 - 1/2
            eta =  2.5 * (psi + lambda);
                                  
            % RUN MODEL

            Yaws_Model_File_baseV2

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
                       
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = log(binopdf(round(13490*ActivePrevvec(i)),13490,0.024));
            LactivePrev_new = log(binopdf(round(13490*ActivePrevvec(i+1)),13490,0.024));

            LlatentPrev_old = log(betapdf((ActivePrevvec(i)/LatentPrevvec(i)),9.78,19.56));
            LlatentPrev_new = log(betapdf((ActivePrevvec(i+1)/LatentPrevvec(i+1)),9.78,19.56));
          
            Ltot_old = LactivePrev_old+LlatentPrev_old;
            Ltot_new = LactivePrev_new+LlatentPrev_new;

            if (Ltot_new - Ltot_old) >= 0
                move_betaONE(i) = 1;
                   
                    Ltot(i+1) = Ltot_new;                  
            else
                y = rand;   

                if y < exp(Ltot_new - Ltot_old)
                    move_betaONE(i) = 1;
                                       
                    Ltot(i+1) = Ltot_new; 
                    
                else
                    betaONEvect(i+1) = betaONEvect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_betaONE(i) = 0;
                                        
                    Ltot(i+1) = Ltot_old;
                end
            end
        end
    end
    
   %% Epsilon 
     if find(Bb == i) > 0   

        triang2 = makedist('Triangular','a',epsilonvect(i)-step_epsilon,'b',epsilonvect(i),'c',epsilonvect(i)+step_epsilon);
        %rng('default'); % For reproducibility
        epsilonvect(i+1) = random(triang2,1,1);
        
        %epsilonvect(i+1) = epsilonvect(i) - step_epsilon + rand(1)*2*step_epsilon;
        
        betaONEvect(i+1)    = betaONEvect(i);
        
            epsilon = epsilonvect(i+1);
            psi =  (3.5/((33 - (0.905/epsilon))*35));
            lambda =  9 * psi;         %1/20 - 1/2
            eta =  2.5 * (psi + lambda);
        
            if epsilonvect(i+1) < 0 || (33 - (0.905/epsilonvect(i+1)))*35 <=0 || eta < 1/60 || epsilonvect(i+1)>20; 
            
            epsilonvect(i+1) = epsilonvect(i);
            move_epsilon(i) = 0;
           
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            LactivePrev(i+1) = LactivePrev_old;
            LlatentPrev(i+1) = LlatentPrev_old;
            Ltot(i+1) = Ltot_old;
           
            else %this is if it falls within the acceptable range
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            epsilon = epsilonvect(i+1);
            psi =  (3.5/((33 - (0.905/epsilon))*35));
            lambda =  9 * psi;         %1/20 - 1/2
            eta =  2.5 * (psi + lambda);
          
            % RUN MODEL

            Yaws_Model_File_baseV2

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
           
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = log(binopdf(round(13490*ActivePrevvec(i)),13490,0.024));
            LactivePrev_new = log(binopdf(round(13490*ActivePrevvec(i+1)),13490,0.024));

            LlatentPrev_old = log(betapdf((ActivePrevvec(i)/LatentPrevvec(i)),9.78,19.56));
            LlatentPrev_new = log(betapdf((ActivePrevvec(i+1)/LatentPrevvec(i+1)),9.78,19.56));
          
            Ltot_old = LactivePrev_old+LlatentPrev_old;
            Ltot_new = LactivePrev_new+LlatentPrev_new;
            
                if (Ltot_new - Ltot_old) >= 0
                    move_epsilon(i) = 1;

                        Ltot(i+1) = Ltot_new;                  
                else
                    y = rand;   

                    if y < exp(Ltot_new - Ltot_old)
                        move_epsilon(i) = 1;

                        Ltot(i+1) = Ltot_new; 
                    else
                    epsilonvect(i+1) = epsilonvect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                                    
                    Ltot(i+1) = Ltot_old;
                    move_epsilon(i) = 0;                   
                    end
                end   
            end
     end
end

%% identify best fit parameters
%ind = find(LactivePrev == max(LactivePrev));
ind = find(Ltot == max(Ltot));

bestfitBetaONE = betaONEvect(ind);
bestfitepsilon = epsilonvect(ind);

%% plot parameters
subplot(3,2,1),plot(betaONEvect),title('Beta')
subplot(3,2,2),plot(epsilonvect),title('Epsilon')
subplot(3,2,4),plot(ActivePrevvec),title('Active Prev')
subplot(3,2,5),plot(LatentPrevvec),title('Latent Prev')
subplot(3,2,6),plot(Ltot),title('Likelihood')
mean(move_betaONE)*2
mean(move_epsilon)*2
%save mcmc_results_v3.mat
toc