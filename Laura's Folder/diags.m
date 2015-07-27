%% plot parameters
subplot(3,2,1),plot(betaONEvect),title('Beta')
subplot(3,2,2),plot(epsilonvect),title('Epsilon')
subplot(3,2,3),plot(gammavect),title('gamma')
subplot(3,2,4),plot(etavect),title('eta')
subplot(3,2,5),plot(thetavect),title('theta')
subplot(3,2,6),plot(lambdavect),title('lambda')

%% get acceptance ratios

mean(move_betaONE)*6
mean(move_epsilon)*6
mean(move_gamma)*6
mean(move_eta)*6
mean(move_theta)*6
mean(move_lambda)*6