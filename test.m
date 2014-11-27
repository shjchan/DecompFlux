%FD in test.mat includes the flux distribution to be decomposed and the
%network.
%EFMdecomp  in test.mat is a decomposition solved using cplex 11 by the
%author.

load test

[EFM1,~,~,~,~,~] = decompflux(FD,1,1);
%CT=1, biomass reaction as the objective; SM=1, suppress calculation message
w1=(EFM1\FD.flux);%weights of the EFMs
disp('Results using CT=1 (EFM1)');
disp(' ');
disp(['Weight of each EFM (w1):   ' num2str(w1')])
disp(' ');
disp(['Absolute sum of EFM1*w1-flux:   ' num2str(sum(abs(EFM1*w1-FD.flux)))]);
%the error of the decomposition

[EFM2,~,~,~,~,~] = decompflux(FD,2,1);
%CT=2, no objective; SM=1, suppress calculation message
w2=(EFM2\FD.flux);%weights of the EFMs
disp('Results using CT=2 (EFM2)');
disp(' ');
disp(['Weight of each EFM (w2):   ' num2str(w2')])
disp(' ');
disp(['Absolute sum of EFM2*w2-flux:   ' num2str(sum(abs(EFM2*w2-FD.flux)))]);
%the error of the decomposition

[~,ind]=min(FD.obj);
Con1=EFM1(ind,:).*w1';%contribution of EFM1 to biomass reaction;
Con2=EFM2(ind,:).*w2';%contribution of EFM2 to biomass reaction;

%Draw histograms to compare the distribution of the contributions of EFMs
%to biomass reaction
xbin=(min([Con1 Con2]):(max([Con1 Con2])-min([Con1 Con2]))/10:max([Con1 Con2]));
x1=min([Con1 Con2])-((max([Con1 Con2])-min([Con1 Con2]))/10);
x2=max([Con1 Con2])+((max([Con1 Con2])-min([Con1 Con2]))/10);
f1=hist(Con1,xbin);
f2=hist(Con2,xbin);
hist(Con1,xbin);
xlabel('Contribution to biomass of the set EFM1 found by CT=1');
ylabel('No. of EFMs');
axis([x1,x2,0,max([f1 f2])]);
figure
hist(Con2,xbin);
xlabel('Contribution to biomass of the set EFM2 found by CT=2');
ylabel('No. of EFMs');
axis([x1,x2,0,max([f1 f2])]);
clear f1 f2 ind x1 x2 xbin


