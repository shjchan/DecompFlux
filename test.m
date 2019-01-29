% flux in example.mat is the flux distribution to be decomposed 
% CbModel in example.mat is the COBRA model (E. coli iAF1260)
% EFMdecomp  in example.mat is a decomposition solved using cplex 11

load('example.mat')

solver = 'cobra';

% CT=1, biomass reaction as the objective; SM=1, suppress calculation message
EFM1 = decompFlux(CbModel, flux, struct('CT', 1, 'solver', solver));
w1 = (EFM1 \ flux);%weights of the EFMs
fprintf('Error of decomposition 1:  %.4e\n', sum(abs(EFM1 * w1 - flux)));

% CT=2, no objective; SM=1, suppress calculation message
EFM2 = decompFlux(CbModel, flux, struct('CT', 2, 'solver', solver));
w2 = (EFM2 \ flux);%weights of the EFMs
fprintf('Error of decomposition 2:  %.4e\n', sum(abs(EFM2*w2-flux)));

% the index for the biomass reaction
ind = find(CbModel.c);
Con1 = EFM1(ind,:) .* w1';  % contribution of EFM1 to biomass reaction;
Con2 = EFM2(ind,:) .* w2';  % contribution of EFM2 to biomass reaction;

%Draw histograms to compare the distribution of the contributions of EFMs
%to biomass reaction
xbin = (min([Con1 Con2]):((max([Con1 Con2]) - min([Con1 Con2]))/10):max([Con1 Con2]));
x1 = min([Con1 Con2])-((max([Con1 Con2])-min([Con1 Con2]))/10);
x2 = max([Con1 Con2])+((max([Con1 Con2])-min([Con1 Con2]))/10);
subplot(2, 1, 1)
f1 = hist(Con1, xbin);
f2 = hist(Con2, xbin);
hist(Con1,xbin);
xlabel('Contribution to biomass of the set EFM1 found by CT=1');
ylabel('No. of EFMs');
axis([x1,x2,0,max([f1 f2])]);
subplot(2, 1, 2)
hist(Con2,xbin);
xlabel('Contribution to biomass of the set EFM2 found by CT=2');
ylabel('No. of EFMs');
axis([x1,x2,0,max([f1 f2])]);
clear f1 f2 ind x1 x2 xbin