function [yn,p1,p2] = isDC(fm,ST,lb0,ub0,flux,k,solver,param)
%[yn,p1,p2] = DC2(fm,ST,lb0,ub0,flux)
%Check the decomposability of flux mode fm 
%by constraining each reaction to have zero flux respectively.
%
%Input:
%fm: flux mode
% ST: stoichiometric matrix
% lb0: lower bound
% ub0: upper bound
% flux: original flux distribution
%solver (optional): solver used to solve LP: 'cplex', 'cobra' or 'matlab'
%                   (extremely slow when using matlab solver)
%                   search automatically if not provided.
%param (optional): optimization parameters
%
%output
%yn: EFM or not
%p1,p2: decomposed flux modes scale to flux
eps0=10^(1); %required sum of fluxes so that flux distribution of all zeros is excluded.
eps7=10^(-30);
eps8=10^(-6);
eps9=10^(-4);
%epsB=10^(-8);
%cor=[];
nz=(fm>eps8);
[Nm,Nr]=size(ST);
indr=(1:Nr)';
indnz=indr(nz);
ubp=ub0;
ubp(fm<=eps8)=0;
if strcmp(solver,'cplex') || strcmp(solver,'matlab')
    f=zeros(Nr,1);
    Aineq=[ST;-ST;-ones(1,Nr)];
    bineq=[eps7*ones(2*Nm,1);-eps0];
end
if strcmp(solver,'cobra')
    LPproblem.A=[ST;-ST;-ones(1,Nr)];
    LPproblem.b=[eps7*ones(2*Nm,1);-eps0];
    LPproblem.c=zeros(Nr,1);
    LPproblem.lb=lb0;
    LPproblem.osense=1;
    LPproblem.csense=char(ones(1,numel(LPproblem.b))*'L');
end
if strcmp(solver,'cplex')
    if ~exist('param','var')
        param=cplexoptimset;
    end
    param.TolFun=10^(-8);
    param.TolRLPFun=10^(-8);
    param.TolXInteger=10^(-8);
end
yn=true;
ct=0;
%disp({'No. of nonzeros:' sum(nz)});
while yn && ct<sum(nz)
    ct=ct+1;
    ubpc=ubp;
    ubpc(indnz(ct))=0;
    if strcmp(solver,'cplex')
        [p1,~,~,out]=cplexlp(f,Aineq,bineq,[],[],lb0,ubpc,[],param);
        if out.cplexstatus==1
            yn=false; %expect no feasible solution. If feasible, then not EFM
        end
    elseif strcmp(solver,'cobra')
        LPproblem.ub=ubpc;
        if exist('param','var')
            sol=solveCobraLP(LPproblem,param);
        else
            sol=solveCobraLP(LPproblem);
        end
        p1=sol.full;
        if sol.stat==1
            yn=false;
        end
    elseif strcmp(solver,'matlab')
        param.Display='off';
        param.Algorithm='simplex';
        [p1,~,exitflag] = linprog(f,Aineq,bineq,[],[],lb0,ubpc,[],param);
        if exitflag==1
            yn=false; %expect no feasible solution. If feasible, then not EFM
        end
    else
        error('A supported solver (Cplex, Cobra or Matlab linprog solver) cannot be identified.')
    end
    %disp(ct);
end
if yn        
    %disp('Non-decomposable. EFM');
    p1=[];
    p2=[];
else
    nz=(p1>eps8);
    p2=fm-min(fm(nz)./p1(nz))*p1;
    p1=min(flux(nz)./p1(nz))*p1;
    p2(p2<=eps9)=0;
    nz2=(p2>0);
    p2=min(flux(nz2)./p2(nz2))*p2;
    imb=sum(abs(ST*p1));
    if imb>eps9
        yn=true;
    %else
        %disp({'imbalance of p1:' imb});
        %if imb>epsB
         %   cor=p1;
         %   [p1] = correct(ST,p1,lb0,flux);
        %end
       % imb=sum(abs(ST*p2));
        %disp({'imbalance of p2:' imb});
       % if imb>epsB
       %     cor=[cor p2];
       %     [p2] = correct(ST,p2,lb0,flux);
       % end
    else
        disp(['The ' num2str(k) '-th EFM is decomposable, not an EFM']);
    end
end

   


end

