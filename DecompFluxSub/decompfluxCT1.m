function [EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT1...
    (ST,flux0,big,lb0,ub0,coeff,block,SM,solver,param)
%[EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT1(ST,flux0,big,lb0,ub0,coeff,block)
%Decompose a flux distribution into a set of EFM
%In the subproblem DC, "p_j<=M(v_j)(a_j)" is used where M is a big number v_j is
%the flux value of the vector 'flux'
%
%Input:
%ST:    m by n stoichiometric matrix 
%       (every reversible reactions should have divided into two irreversible reactions in the matrix)
%flux:  n by 1 flux distribution vector
%big:   the large number M used in the optimization problem 
%       (may need to adjust so that it's not too low or high, 
%       because when it is too large, M(v_j) becomes too large for some
%       reaction and Cplex may not tolerate. When too low, M(v_j) may
%       becomes not large enough for some reactions and some real solutions
%       may be ruled out. 10^8~10^9 should be good enough)
%lb0:   n by 1 lower bound vector (default zeros if empty)
%ub0:   n by 1 upper bound vector (for scaling purpose, default the large number M if empty)
%coeff: n by 1 objective coefficient vector
%block: If you have k flux modes to be excluded from the possible solutions,
%       then you have a n by k matrix in which each column has 1 for reactions 
%       with non-zero fluxes in the flux mode and 0 for other entries.
% eg. in a network with three reactions only,
%     if block =
%           [1 0;...
%            1 1;...
%            0 1];  
%     then the flux mode with non-zero fluxes in reaction 1 and 2, and the
%     flux mode with non-zero fluxes in reaction 2 and 3, are being blocked
%     from being the solution of the optimization problem.
%SM: level of display, 0 (no display) - 3 (most detailed)
%solver: currently supported solver: 'cobra' or 'cplex'
%param: optimization parameter, for 'cobra' only
%
%Output:
%EFM: the set of EFM
%FM: the set of flux modes found during computation
%Rshut: reactions being shut down when computing the n-th flux mode
%Reduce: the sum of fluxes (of the flux distribution to be decomposed) that
%       is subtracted after finding the n-th EFM
%Corr:  2 by N cell, if corrections of flux modes have taken place due to 
%       high level of flux imbalance, the no. of flux modes and the flux
%       modes before correction will be recorded in Corr.
%depth: When an EFM is found, the depth of the stack, ie., 'current node'
%       number and the depth of the stack after removing useless flux modes, 
%       will be recorded in this N by 2 matrix 'depth'.
%
%Internal parameters: (need to edit this script if you want to change)
eps0=max(flux0)/10^(7); %tolerance for zero in processing the flux modes in the stack, 
%ie. entries<eps0 will be treated as strict zeros.
epsPW=10^(-10);%tolerance for zero in processing EFM found.
epsPW2=10^(-7);%tolerance for zero after the found EFM is scaled up to the flux distribution


ts=5000;
tsefm=ts;
tsfm=ts;
tscur=ts;
nzf=(flux0>eps0);
flux=flux0*min(ub0(nzf)./flux0(nzf));
art=flux;
Nr=length(flux);
if isempty(lb0)
    lb0=zeros(Nr,1);
end
if isempty(ub0)
    ub0=big*ones(Nr,1);
end
PWtree=zeros(Nr,ts);
EFM=zeros(Nr,ts);
Rshut=cell(ts,1);
Reduce=zeros(ts,1);
FM=zeros(Nr,ts);
FM(:,1)=art;
depth=zeros(ts,2);
bt=0;
tsbt=ts;
cur=1;
PWtree(:,cur)=art;
Nefm=0;
Nfm=1;
crmax=100;
Ncr=100;
Corr=cell(2,crmax);
cr=0;
reduceP = 0;
fluxP = flux0;
while cur>0
    if SM >= 2
        disp(['Current node:' num2str(cur) '. No. of nodes passed:' num2str(Nfm)]);
    end
    if cur==1 && min(art(art~=0))<0.1
        ratio=max([min(art(art~=0));max(art./ub0)]);
        art=art/ratio;
    end
    if sum(art>eps0)>0
        if strcmp(solver,'cplex')
            [yn,p2,Rz,Cor] = decompfluxCT1sub(ST,art,lb0,PWtree(:,1),big,coeff,block,SM);
        elseif strcmp(solver,'cobra')
            if exist('param','var')
                [yn,p2,Rz,Cor] = decompfluxCT1subCOBRA(ST,art,lb0,PWtree(:,1),big,coeff,block,SM,param);
            else
                [yn,p2,Rz,Cor] = decompfluxCT1subCOBRA(ST,art,lb0,PWtree(:,1),big,coeff,block,SM);
            end
        else
            error('A supported solver (Cplex or Cobra solver) cannot be identified.')
        end
        if ~isempty(Cor)
            if cr>=Ncr
                Corr=[Corr cell(2,crmax)];
                Ncr=crmax+Ncr;
            end
            cr=cr+1;
            Corr{1,cr}=Nfm;
            Corr{2,cr}=Cor;
        end
        if yn
            Nefm=Nefm+1;
            if Nefm>tsefm
                EFM=[EFM zeros(Nr,ts)];
                Reduce=[Reduce;zeros(ts,1)];
                tsefm=tsefm+ts;
            end
            EFM(:,Nefm)=art;%PWtree(:,cur);
            efmc=PWtree(:,cur);
            if cur==1
                if SM ~= 0
                    disp('Algorithm is terminated.');
                end
                break;
            else
                in=diag(1./efmc(efmc>0));
                ra=in*PWtree(efmc>0,1:cur-1);
                ra=min(ra,[],1);
                subtract=efmc*ra;
                PWtree(:,1:cur-1)=PWtree(:,1:cur-1)-subtract;
                PWtree(PWtree<=epsPW)=0;
                Reduce(Nefm)=sum(subtract(:,1));
                PWbin=(PWtree(:,1:cur-1)>epsPW);
                PWflux=logical((PWtree(:,1)<=epsPW)*ones(1,cur-1));
                res=(sum(PWbin&PWflux)==0);
                bt=bt+1;
                if bt>tsbt
                    depth=[depth;zeros(ts,2)];
                    tsbt=tsbt+ts;
                end
                depth(bt,:)=[cur sum(res)];
                cur=sum(res);
                PWtree(:,1:cur)=PWtree(:,res);
                ratio2=max(diag(1./ub0)*PWtree(:,1:cur));
                ratio2=1./ratio2;
                PWtree(:,1:cur)=PWtree(:,1:cur)*diag(ratio2);
                PWtree(PWtree<=epsPW2)=0;
                art=PWtree(:,cur);
                while sum(art>eps0)==0
                    cur=cur-1;
                    if cur==0
                        break
                    end
                    art=PWtree(:,cur);
                end
                if SM~=0
                    r = min(fluxP(efmc~=0)./efmc(efmc(:)~=0));
                    fluxP = fluxP - efmc*r;
                    fluxP(fluxP<epsPW)=0;
                    reduceP = reduceP + sum(efmc)*r*100/sum(flux0);
                    clockN = clock;
                    if clockN(5) < 10
                        c5 = ['0' num2str(clockN(5))];
                    else
                        c5 = num2str(clockN(5));
                    end
                    disp(['No. of EFMs found: ' num2str(Nefm) ...
                        '. Flux decomposed: ' num2str(round(reduceP)) ...
                        ' %. ' date ' , ' num2str(clockN(4)) ':' c5]);
                end
            end
            if cur==0
                if SM ~= 0
                    disp('Algorithm is terminated.');
                end
                break;
            end
            Nfm=Nfm+1;
            FM(:,Nfm)=art;
        else%branching
            if cur+1>tscur
                PWtree=[PWtree zeros(Nr,ts)];
                tscur=tscur+ts;
            end
            Nfm=Nfm+1;
            if Nfm>tsfm
                FM=[FM zeros(Nr,ts)];
                Rshut=[Rshut;cell(ts,1)];
                tsfm=tsfm+ts;
            end
            Rshut{Nfm}=Rz;
            p2(p2<0)=0;%rectify p2
            FM(:,Nfm)=p2;
            PWtree(:,cur+1)=p2;
            cur=cur+1;
            art=p2;
        end
    else
        if SM ~= 0
            disp('Looping');
        end
        if cur==1
            if SM ~= 0
                disp('Algorithm is terminated.');
            end
            Nfm=Nfm+1;
            FM(:,Nfm)=art;
            break;
        else
            while sum(art>eps0)==0 && cur>1
                cur=cur-1;
                art=PWtree(:,cur);
                Nfm=Nfm+1;
                FM(:,Nfm)=art;
            end
        end
        if cur==0
            if SM ~= 0
                disp('error.');
            end
            break;
        end
    end
    clear p2;
end
EFM=EFM(:,1:Nefm);
w = EFM\flux0;
EFM = EFM * diag(w);
FM=FM(:,1:Nfm);
for j2 = 1:size(FM,2);
    r = min(flux0(FM(:,j2)~=0)./FM(FM(:,j2)~=0,j2));
    FM(:,j2) = FM(:,j2) * r;
end
Reduce=Reduce(1:Nefm) .* w;
Rshut=Rshut(1:Nfm,:);
Corr=Corr(:,1:cr);
depth=depth(1:bt,:);
%disp({'No. of EFM:' Nefm;'Flux by reduction:' sum(Reduce)});
end

