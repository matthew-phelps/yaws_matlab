%plots

subplot(4,2,1),plot(betaONEvect),title('Beta Vect')
subplot(4,2,2),plot(epsilonvect),title('Epsilon Vect')
subplot(4,2,3),plot(IP1./N),title('IP')
subplot(4,2,4),plot(L1./N),title('L1')
subplot(4,2,5),plot(IS./N),title('IS')
subplot(4,2,6),plot(L2./N),title('L2')
subplot(4,2,7),plot(IT./N),title('IT')
subplot(4,2,8),plot(R./N),title('R')