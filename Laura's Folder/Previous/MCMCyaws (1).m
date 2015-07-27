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

        
loops2 = 200; %200150;

ActivePrevvec = zeros(loops2+1,1); %creating space to be populated with the prevelance values from the model. Loops is in reference to how many iterations we will be doing

LactivePrev = zeros(loops2+1,1); %This is space where likelihood values will be stored

YawsParameterFile_aug % my parameter file

Yaws_Model_File %my model file

ActivePrevvec(1,1) = ActivePrev(end); %we only save values that are accepted by the Hastings algorithm. 

%% These are the three parameters we want to fit

betaONEvect = zeros(loops2+1,1);
betaONEvect(1) = betaONE;      % 0.000001 - 0.00001
move_betaONE = zeros(loops2,1);
step_betaONE = 0.01*betaONEvect(1);

gammavect = zeros(loops2+1,1); %this is the blank matrix
gammavect(1) = gamma;            %1/9 - 1/3                                      % point estimate that we start with. we generate a vector of accepted values for parameter - this is the starting one
move_gamma = zeros(loops2,1); %setting up space for dummy variable - if it moves it gets a one, if not, a zero
step_gamma = 0.1*gammavect(1); %step size. Can make it dependent on initial value

epsilonvect = zeros(loops2+1,1);
epsilonvect(1) = epsilon;          %1/24 - 0
move_epsilon = zeros(loops2,1);
step_epsilon = 0.01*epsilonvect(1);

etavect = zeros(loops2+1,1);
etavect(1) = eta;             %1/60-1/3
move_eta = zeros(loops2,1);
step_eta = 0.01*etavect(1); 

thetavect = zeros(loops2+1,1);
thetavect(1) = theta;           %1/6 - 1/0.5
move_theta = zeros(loops2,1);
step_theta = 0.01*etavect(1);

lambdavect = zeros(loops2+1,1);
lambdavect(1) = lambda;           %1/20 - 1/2
move_lambda = zeros(loops2,1);
step_lambda = 0.05*etavect(1);


%% LIKELIHOOD FUNCTIONS

LactivePrev_old = binopdf(round(1176*ActivePrevvec(1)),1176,0.047);

% the probability (P) comes from the data, the x comes from the model, N comes
% from the study population from which P was estimated
% N = 1176; P=.047

%Take the x output from the model and use that as the new input
%(round(1176*prev), 1776, 0.047)

LactivePrev(1) = LactivePrev_old;

Aa = [1:6:200150];%the big number at the end is the total number of loops we do
Bb = [2:6:200150];
Cc = [3:6:200150];
Dd = [4:6:200150];
Ee = [5:6:200150];
Ff = [6:6:200150];
%% For loops

for i = 1:loops2 
    if find(Aa == i) > 0    %sequence, administrative thing, making sure that only 1 variable is moved at a time
                            %If i=1 we go through this first loop, if i=2
                            %we go through the second loop
                            
        %taking intial value of param and adding a random # (drawn from a uniform dist between 0 and 1)
        %to it. The # can be pos or neg, but does not exceed step size specified above
        gammavect(i+1) = gammavect(i) - step_gamma + rand(1)*2*step_gamma;
        
        %Here we specify that we are leaving all other parameters the same
        %value they had in the previous step
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
        
        %here we evaluate the new gamma param in the context of the
        %existing values of the other parameters
        
        %if the random value falls outside of the range below we will not use
        %it, we only use values within the range below to be considered.
        %Since we know something about the duration of IP (b/w 3 and 6 
        % months), we can use that as a starting point for this range. 
        if gammavect(i+1) <= 1/9 || gammavect(i+1) > 1/3; %

            % if it falls in rejection region specified above we go back to
            % the previous step value for gamma
            gammavect(i+1) = gammavect(i);
            
            %for the dummy variable (move_gamma) we assigns a 0, which
            %means no movement
            move_gamma(i) = 0;
            
            %if new gamma is rejected we revert to the previous value for
            %ActivePrevvec
            ActivePrevvec(i+1) = ActivePrevvec(i);

            LactivePrev(i+1) = LactivePrev_old;
           
        else %this is if it falls within the acceptable range
            betaONE = betaONEvect(i+1);
            betaTWO = betaONE;
            epsilon = epsilonvect(i+1);
            eta = etavect(i+1);
            theta = thetavect(i+1);
            lambda = lambdavect(i+1);
            gamma = gammavect(i+1);
            
          
            % RUN MODEL

            %yaws model file
            Yaws_Model_File

            ActivePrevvec(i+1) = ActivePrev(end);
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_gamma(i) = 1;
                   
                    % store the new likelihood value if the new likelihood
                    % is greater than the old likelihood
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;   %even if new L is less than old, we still want
                            %to accept the new value a certain random percentage of times

                if y < (LactivePrev_new / LactivePrev_old)
                    move_gamma(i) = 1;
                    
                    
                    LactivePrev(i+1) = LactivePrev_new;
                   

                else
                    gammavect(i+1) = gammavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_gamma(i) = 0;
                    
                    
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
    end
    if find(Bb == i) > 0    %sequence, administrative thing, making sure that only 1 variable is moved at a time
                            %If i=1 we go through this first loop, if i=2
                            %we go through the second loop
                            
        %taking intial value of param and adding a random # (drawn from a uniform dist between 0 and 1)
        %to it. The # can be pos or neg, but does not exceed step size specified above
        epsilonvect(i+1) = epsilonvect(i) - step_epsilon + rand(1)*2*step_epsilon;
        
        %Here we specify that we are leaving all other parameters the same
        %value they had in the previous step
        gammavect(i+1)      = gammavect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
        
        %here we evaluate the new gamma param in the context of the
        %existing values of the other parameters
        
        %if the random value falls outside of the range below we will not use
        %it, we only use values within the range below to be considered.
  
        if epsilonvect(i+1) <= 0 || epsilonvect(i+1) > 1/24; 

            % if it falls in rejection region specified above we go back to
            % the previous step value for gamma
            epsilonvect(i+1) = epsilonvect(i);
            
            %for the dummy variable (move_gamma) we assigns a 0, which
            %means no movement
            move_epsilon(i) = 0;
            
            %Does this stay the same when doing this variable??
            ActivePrevvec(i+1) = ActivePrevvec(i);

            LactivePrev(i+1) = LactivePrev_old;
           
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
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_epsilon(i) = 1;
                   
                    % store the new likelihood value if the new likelihood
                    % is greater than the old likelihood
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;   %even if new L is less than old, we still want
                            %to accept the new value a certain random percentage of times

                if y < (LactivePrev_new / LactivePrev_old)
                    move_epsilon(i) = 1;
                    
                    
                    LactivePrev(i+1) = LactivePrev_new;
                   

                else
                    epsilonvect(i+1) = epsilonvect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_epsilon(i) = 0;
                    
                    
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
    end
   if find(Cc == i) > 0    
                            
        etavect(i+1) = etavect(i) - step_eta + rand(1)*2*step_eta;
        
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        betaONEvect(i+1)    = betaONEvect(i);
        
      
        if etavect(i+1) <= 1/60 || etavect(i+1) > 1/3; %Ask Laura about range for eta

            etavect(i+1) = etavect(i);
            
            move_eta(i) = 0;
            
            ActivePrevvec(i+1) = ActivePrevvec(i);

            LactivePrev(i+1) = LactivePrev_old;
           
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
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_eta(i) = 1;
                   
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;   

                if y < (LactivePrev_new / LactivePrev_old)
                    move_eta(i) = 1;
                                        
                    LactivePrev(i+1) = LactivePrev_new;
                   
                else
                    etavect(i+1) = etavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_eta(i) = 0;
                                        
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
   end

    %############################
    %############################
    %  Theta  ###################
    
    if find(Dd == i) > 0    
           
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

            LactivePrev(i+1) = LactivePrev_old;
           
        else 
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
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_theta(i) = 1;
                   
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;  

                if y < (LactivePrev_new / LactivePrev_old)
                    move_theta(i) = 1;
                    
                    
                    LactivePrev(i+1) = LactivePrev_new;
                   

                else
                    thetavect(i+1) = thetavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_theta(i) = 0;
                    
                    
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
    end
    
    %#############################
    %#############################
    %  Lambda  ###################
    
    if find(Ee == i) > 0    
               
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

            LactivePrev(i+1) = LactivePrev_old;
           
        else 
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
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_lambda(i) = 1;
                   
                    % store the new likelihood value if the new likelihood
                    % is greater than the old likelihood
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;   %even if new L is less than old, we still want
                            %to accept the new value a certain random percentage of times

                if y < (LactivePrev_new / LactivePrev_old)
                    move_lambda(i) = 1;
                    
                    
                    LactivePrev(i+1) = LactivePrev_new;
                   

                else
                    lambdavect(i+1) = lambdavect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_lambda(i) = 0;
                    
                    
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
    end
   
       %#############################
    %#############################
    %  betaONE  ###################
    
    if find(Ff == i) > 0    
                            
        %taking intial value of param and adding a random # (drawn from a uniform dist between 0 and 1)
        %to it.
        betaONEvect(i+1) = betaONEvect(i) - step_betaONE + rand(1)*2*step_betaONE;
        
        %Here we specify that we are leaving all other parameters the same
        %value they had in the previous step
        gammavect(i+1)      = gammavect(i);
        epsilonvect(i+1)    = epsilonvect(i);
        etavect(i+1)        = etavect(i);
        thetavect(i+1)      = thetavect(i);
        lambdavect(i+1)     = lambdavect(i);
        
      
        if betaONEvect(i+1) <= .0002 || betaONEvect(i+1) > .002; %Ask Laura about range for epsilon

            % if it falls in rejection region specified above we go back to
            % the previous step value for gamma
            betaONEvect(i+1) = betaONEvect(i);
            
            move_betaONE(i) = 0;
            
            %Does this stay the same when doing this variable??
            ActivePrevvec(i+1) = ActivePrevvec(i);

            LactivePrev(i+1) = LactivePrev_old;
           
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
           
            
            % LIKELIHOOD FUNCTIONS

            LactivePrev_old = binopdf(round(1176*ActivePrevvec(i)),1176,0.047);

            LactivePrev_new = binopdf(round(1176*ActivePrevvec(i+1)),1176,0.047);

            % here we compare the two likelihoods
            if (LactivePrev_new / LactivePrev_old) >= 1
                move_betaONE(i) = 1;
                   
                    % store the new likelihood value if the new likelihood
                    % is greater than the old likelihood
                    LactivePrev(i+1) = LactivePrev_new;
                  
            % if new L is less than the old
            else
                y = rand;   %even if new L is less than old, we still want
                            %to accept the new value a certain random percentage of times

                if y < (LactivePrev_new / LactivePrev_old)
                    move_betaONE(i) = 1;
                    
                    
                    LactivePrev(i+1) = LactivePrev_new;
                   

                else
                    betaONEvect(i+1) = betaONEvect(i);
                    ActivePrevvec(i+1) = ActivePrevvec(i);
                    
                    move_betaONE(i) = 0;
                    
                    
                    LactivePrev(i+1) = LactivePrev_old;
                   
                end
            end   
        end
    end 
end

%% identify best fit parameters
ind = find(LactivePrev == max(LactivePrev));

bestfitEpsilon = epsilonvect(ind);
bestfitEta = etavect(ind);
bestfitGamma = gammavect(ind);
bestfitTheta = thetavect(ind);
bestfitLambda = lambdavect(ind);
bestfitBetaONE = betaONEvect(ind);



cend = fix(clock);