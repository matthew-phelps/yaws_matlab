%% plot parameters
subplot(3,2,1),plot(betaONEvect),title('Beta')
subplot(3,2,2),plot(thetavect),title('Theta')
subplot(3,2,3),plot(etavect),title('Eta')
subplot(3,2,4),plot(ActivePrevvec),title('Active Prev')
subplot(3,2,5),plot(LatentPrevvec),title('Latent Prev')
subplot(3,2,6),plot(Ltot),title('Likelihood')

%% get acceptance ratios

mean(move_betaONE)*3
mean(move_theta)*3
mean(move_eta)*3
