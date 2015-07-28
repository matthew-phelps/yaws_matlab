%% plot parameters
subplot(3,2,1),plot(betaONEvect),title('Beta')
subplot(3,2,2),plot(epsilonvect),title('Epsilon')
subplot(3,2,4),plot(ActivePrevvec),title('Active Prev')
subplot(3,2,5),plot(LatentPrevvec),title('Latent Prev')
subplot(3,2,6),plot(Ltot),title('Likelihood')

%% get acceptance ratios

mean(move_betaONE)*2
mean(move_epsilon)*2

