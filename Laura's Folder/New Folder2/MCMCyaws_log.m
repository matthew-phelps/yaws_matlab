cbegin = fix(clock); %this gives you the time it started and time it ended

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
        
loops2 = 2000; %200150;

ActivePrevvec = zeros(loops2+1,1); %creating space to be populated with the prevalance values from the model. Loops is in reference to how many iterations we will be doing

LatentPrevvec = zeros(loops2+1,1); % remove if we don't use ratio method 

LactivePrev = zeros(loops2+1,1); %This is space where likelihood values will be stored
LlatentPrev = zeros(loops2+1,1);

Ltot = zeros(loops2+1,1); % this is the likelihood of both active and latent combined

YawsParameterFile_aug % my parameter file

Yaws_Model_File %my model file

ActivePrevvec(1) = ActivePrev(end); %we only save values that are accepted by the Hastings algorithm. 

LatentPrevvec(1) = LatentPrev(end); % remove if we don't use ratio method

%% These are the three parameters we want to fit

betaONEvect = zeros(loops2+1,1);
betaONEvect(1) = betaONE;      % 0.000001 - 0.00001
move_betaONE = zeros(loops2,1);
step_betaONE = 0.2*betaONEvect(1);

gammavect = zeros(loops2+1,1); %this is the blank matrix
gammavect(1) = gamma;            %1/9 - 1/3                                      % point estimate that we start with. we generate a vector of accepted values for parameter - this is the starting one
move_gamma = zeros(loops2,1); %setting up space for dummy variable - if it moves it gets a one, if not, a zero
step_gamma = 0.1*gammavect(1); %step size. Can make it dependent on initial value

epsilonvect = zeros(loops2+1,1);
epsilonvect(1) = epsilon;          %1/24 - 0
move_epsilon = zeros(loops2,1);
step_epsilon = 0.2*epsilonvect(1);

etavect = zeros(loops2+1,1);
etavect(1) = eta;             %1/60-1/3
move_eta = zeros(loops2,1);
step_eta = 0.0001*etavect(1); 

thetavect = zeros(loops2+1,1);
thetavect(1) = theta;           %1/6 - 1/0.5
move_theta = zeros(loops2,1);
step_theta = 0.1*thetavect(1);

lambdavect = zeros(loops2+1,1);
lambdavect(1) = lambda;           %1/20 - 1/2
move_lambda = zeros(loops2,1);
step_lambda = 0.1*lambdavect(1);

%% LIKELIHOOD FUNCTIONS

LactivePrev_old = log(binopdf(round(1176*ActivePrevvec(1)),1176,0.047));
LlatentPrev_old = log(binopdf(round(1176*LatentPrevvec(1)),1176,(2*0.047))); % remove if we don't use ratio method
Ltot_old = LactivePrev_old + LlatentPrev_old;
% the probability (P) comes from the data, the x comes from the model, N comes
% from the study population from which P was estimated
% N = 1176; P=.047

%Take the x output from the model and use that as the new input
%(round(1176*prev), 1776, 0.047)

LactivePrev(1) = LactivePrev_old;
LlatentPrev(1) = LlatentPrev_old;
Ltot(1) = Ltot_old;

Aa = [1:6:200150];%the big number at the end is the total number of loops we do
Bb = [2:6:200150];
Cc = [3:6:200150];
Dd = [4:6:200150];
Ee = [5:6:200150];
Ff = [6:6:200150];
%% For loops


for i = 1:loops2 
     %% BetaONE
    
    if find(Aa == i) > 0    
                            
        betaONEvect(i+1) = betaONEvect(i) - step_betaONE + rand(1)*2*step_betaONE;
        
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
              
        if betaONEvect(i+1) <= 0 || betaONEvect(i+1) > 20; %0.0002, 0.002

            betaONEvect(i+1) = betaONEvect(i);
            
            move_betaONE(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            LactivePrev(i+1) = LactivePrev_old;
            Ltot(i+1) = Ltot_old;
            
        else 
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            gamma = gammavect(i+1);
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
                                 
            % RUN MODEL

            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
                       
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = log(binopdf(round(1176*ActivePrevvec(i)),1176,0.047));
            LactivePrev_new = log(binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = log((binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047))));
            LlatentPrev_new = log((binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047))));
          
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
    
    %% gamma
    
    if find(Bb == i) > 0                               
        gammavect(i+1) = gammavect(i) - step_gamma + rand(1)*2*step_gamma;
        
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
        
        if gammavect(i+1) <= 1/9 || gammavect(i+1) > 1/.75; %

            gammavect(i+1) = gammavect(i);
            
            move_gamma(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
           
            Ltot(i+1) = Ltot_old;
            
        else %this is if it falls within the acceptable range
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
            gamma = gammavect(i+1);           
          
            % RUN MODEL
            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
                       
            % LIKELIHOOD FUNCTIONS
            LactivePrev_old = (binopdf(round(1176*ActivePrevvec(i)),1176,0.047));
            LactivePrev_new = (binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = (binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047)));
            LlatentPrev_new = (binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047)));
          
            Ltot_old = LactivePrev_old*LlatentPrev_old;
            Ltot_new = LactivePrev_new*LlatentPrev_new;
            
            if (Ltot_new / Ltot_old) >= 1
                move_gamma(i) = 1;
                   
                    Ltot(i+1) = Ltot_new;
                  
            else
                y = rand;  

                if y < (Ltot_new / Ltot_old)
                    move_gamma(i) = 1;
                                       
                    Ltot(i+1) = Ltot_new;
                   
                else
                    gammavect(i+1) = gammavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_gamma(i) = 0;
                                       
                    Ltot(i+1) = Ltot_old;
                   
                end
            end   
        end
    end
    
    %% Epsilon
    
    if find(Cc == i) > 0   
        
        epsilonvect(i+1) = epsilonvect(i) - step_epsilon + rand(1)*2*step_epsilon;
        
        gammavect(i+1)      = gammavect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
  
        if epsilonvect(i+1) < 0 || epsilonvect(i+1) > 1/24; 

            epsilonvect(i+1) = epsilonvect(i);
            
            move_epsilon(i) = 0;
           
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);

            Ltot(i+1) = Ltot_old;
            
        else %this is if it falls within the acceptable range
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            gamma = gammavect(i+1);
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
          
            % RUN MODEL

            %yaws model file
            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
            
            % LIKELIHOOD FUNCTIONS
            LactivePrev_old = (binopdf(round(1176*ActivePrevvec(i)),1176,0.047));
            LactivePrev_new = (binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = (binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047)));
            LlatentPrev_new = (binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047)));
          
            Ltot_old = LactivePrev_old*LlatentPrev_old;
            Ltot_new = LactivePrev_new*LlatentPrev_new;
            
            if (Ltot_new / Ltot_old) >= 1
                move_epsilon(i) = 1;
                   
                    Ltot(i+1) = Ltot_new;
                  
            else
                y = rand;   
                
                if y < (Ltot_new / Ltot_old)
                    move_epsilon(i) = 1;
                    
                    Ltot(i+1) = Ltot_new;

                else
                    epsilonvect(i+1) = epsilonvect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_epsilon(i) = 0;
                                     
                    Ltot(i+1) = Ltot_old;                   
                end
            end   
        end
    end
    
    %% Eta
    
   if find(Dd == i) > 0    
                            
        etavect(i+1) = etavect(i) - step_eta + rand(1)*2*step_eta;
        
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
              
        if etavect(i+1) <= 0 || etavect(i+1) > 3; 
            
            etavect(i+1) = etavect(i);
            
            move_eta(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            Ltot(i+1) = Ltot_old;
           
        else 
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            gamma = gammavect(i+1);
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
            
            % RUN MODEL

            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
                       
            % LIKELIHOOD FUNCTIONS
            LactivePrev_old = log((binopdf(round(1176*ActivePrevvec(i)),1176,0.047)));
            LactivePrev_new = log(binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = log(binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047)));
            LlatentPrev_new = log(binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047)));
          
            Ltot_old = LactivePrev_old+LlatentPrev_old;
            Ltot_new = LactivePrev_new+LlatentPrev_new;
            
            if (Ltot_new - Ltot_old) >= 0
                move_eta(i) = 1;
                   
                    Ltot(i+1) = Ltot_new;
                  
            else
                y = rand;  

                if y < (Ltot_new - Ltot_old)
                    move_eta(i) = 1;
                                       
                    Ltot(i+1) = Ltot_new;
                   
                else
                    etavect(i+1) = etavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_eta(i) = 0;
                                       
                    Ltot(i+1) = Ltot_old;                  
                end
            end   
        end
   end
   
   %% Theta 
   
    if find(Ee == i) > 0    
           
        thetavect(i+1) = thetavect(i) - step_theta + rand(1)*2*step_theta;
        
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
              
        if thetavect(i+1) <= 1/6  || thetavect(i+1) >1/0.5 ; %get more specific numbers

            thetavect(i+1) = thetavect(i);
            
            move_theta(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            LactivePrev(i+1) = LactivePrev_old;
            Ltot(i+1) = Ltot_old;
           
        else 
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            gamma = gammavect(i+1);
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
                      
            % RUN MODEL

            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
           
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = (binopdf(round(1176*ActivePrevvec(i)),1176,0.047));
            LactivePrev_new = (binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = (binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047)));
            LlatentPrev_new = (binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047)));
          
            Ltot_old = LactivePrev_old*LlatentPrev_old;
            Ltot_new = LactivePrev_new*LlatentPrev_new;
            
            if (Ltot_new / Ltot_old) >= 1
                move_theta(i) = 1;
                   
                    Ltot(i+1) = Ltot_new;
            else
                y = rand;  

                if y < (Ltot_new / Ltot_old)
                    move_theta(i) = 1;
                                       
                    Ltot(i+1) = Ltot_new;
                   
                else
                    thetavect(i+1) = thetavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_theta(i) = 0;
                                       
                    Ltot(i+1) = Ltot_old;                   
                end
            end   
        end
    end    
    
    %% Lambda
    
    if find(Ff == i) > 0    
               
        lambdavect(i+1) = lambdavect(i) - step_lambda + rand(1)*2*step_lambda;
        
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
             
        if lambdavect(i+1) <= 1/120 || lambdavect(i+1) > 1/2; 

            lambdavect(i+1) = lambdavect(i);
            
            move_lambda(i) = 0;
                      
            ActivePrevvec(i+1) = ActivePrevvec(i);
            LatentPrevvec(i+1) = LatentPrevvec(i);
            
            Ltot(i+1) = Ltot_old;
            
        else 
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            gamma = gammavect(i+1);
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
                     
            % RUN MODEL

            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
            LatentPrevvec(i+1) = LatentPrev(end);
                      
            % LIKELIHOOD FUNCTIONS
            LactivePrev_old = (binopdf(round(1176*ActivePrevvec(i)),1176,0.047));
            LactivePrev_new = (binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047));

            LlatentPrev_old = (binopdf(round(1176*LatentPrevvec(i)),1176,(2*0.047)));
            LlatentPrev_new = (binopdf(round(1176*LatentPrevvec(i+1)),1176,(2*0.047)));
          
            Ltot_old = LactivePrev_old*LlatentPrev_old;
            Ltot_new = LactivePrev_new*LlatentPrev_new;
            
            if (Ltot_new / Ltot_old) >= 1
                move_lambda(i) = 1;
                   
                Ltot(i+1) = Ltot_new;
                  
            else
                y = rand;   

                if y < (Ltot_new / Ltot_old)
                    move_lambda(i) = 1;
                                      
                    Ltot(i+1) = Ltot_new;
                   
                else
                    lambdavect(i+1) = lambdavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    LatentPrevvec(i+1) = LatentPrevvec(i);
                    
                    move_lambda(i) = 0;
                                        
                    Ltot(i+1) = Ltot_old;                  
                end
            end   
        end
    end
    
    
   
end

%% identify best fit parameters
%ind = find(LactivePrev == max(LactivePrev));
ind = find(Ltot == max(Ltot));

bestfitEpsilon = epsilonvect(ind);
bestfitEta = etavect(ind);
bestfitGamma = gammavect(ind);
bestfitTheta = thetavect(ind);
bestfitLambda = lambdavect(ind);
bestfitBetaONE = betaONEvect(ind);

cend = fix(clock);