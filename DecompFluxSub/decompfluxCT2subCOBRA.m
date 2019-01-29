function [yn,p2,Rz,Cor] = decompfluxCT2subCOBRA(ST,fm,lb0,flux,big,coeff,block,SM,param)
%[yn,p2,Rz,Cor] = decompfluxCT2subCOBRA(ST,fm,lb0,flux,big,coeff)
%Check the decomposability of flux mode fm by MILP using COBRA solver
%Using the constraint "p_j<=M(delta_j)(a_j)" where M is a big number and 
%delta_j is 0 if v_j in 'flux' is zero, otherwise 1.
%
%Input:
%ST:    m by n stoichiometric matrix 
%       (every reversible reactions should have divided into two irreversible reactions in the matrix)
%fm:    n by 1 vector, the flux mode to be checked, only its set of zeros
%       is important, the actual value is not important.
%lb0:   n by 1 lower bound vector
%flux:  n by 1 flux distribution vector, for constraining the problem.
%big:   the large number M used in the optimization problem 
%       (not recommended to be too high or too low for the feasibility of 
%       the optimization problem, 10^6~10^9 should be good enough)
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
%param: optional optimization parameters for COBRA solver
%
%Output:
%yn: true if EFM, false if not
%p2: solution flux mode
%Rz: reaction being shut down
%Cor: solution flux mode before correction if correction has happened.

%Internal parameter:
epsIMB=10^(-30);
%tolerance for flux imbalance in constraints. In our testing with Cplex, actually using
%two inequality constraints with some imbalance tolerance in RHS gives better
%results than using one equality constraint with 0 in RHS.

eps8=max(fm)/10^(7);
%tolerance for zeros when handling fm

eps9=10^(-4);
%tolerance for zeros when handling solutions which is possibly false-positive
%(feasible soultions returned when there should be no solution)

eps10=10^(-7);
%tolerance for imbalance to admit a returned solution is true-positive

epsB=10^(-7);
%tolerance for imbalance to correct a true-positive solution

epsP=10^(-5);
%tolerance for zeros when handling returnd solutions.

[Nm,Nr]=size(ST);
f=[zeros(Nr,1);full(coeff(:))];
%Constraint for blocked EFMs
if isempty(block)
    BK=[];
else
    BK=sparse(zeros(size(block')));
    BK(logical(block'))=1;
    BK=[BK sparse(zeros(size(block,2),Nr))];
end
%Constraint matrix
MILPproblem.A = sparse([[zeros(2*Nm,Nr) [ST;-ST]];...
       [eye(Nr)        -eye(Nr)];...
       [-big*eye(Nr)  eye(Nr)];...
       [ones(1,Nr)  zeros(1,Nr)];... 
       [-ones(1,Nr) zeros(1,Nr)];... 
       BK;...
       ]);
%RHS
if isempty(BK)
    MILPproblem.b=[epsIMB*ones(2*Nm,1);zeros(2*Nr,1);sum(fm>eps8)-1;-1];
else
    MILPproblem.b=[epsIMB*ones(2*Nm,1);zeros(2*Nr,1);sum(fm>eps8)-1;-1;sum(logical(block'),2)-1];
end
MILPproblem.c=f; %objective vector
ubp=big*ones(Nr,1);
ubp(fm<=eps8)=0;
uba=ones(Nr,1);
uba(fm<=eps8)=0;
lba=zeros(Nr,1);
MILPproblem.lb=[lba;lb0(:)]; %lower bound
MILPproblem.ub=[uba;ubp]; %upper bound
MILPproblem.osense=1; %minimization
MILPproblem.csense=char(ones(1,numel(MILPproblem.b))*'L'); %sense of constraints
MILPproblem.vartype=char([ones(1,Nr)*'B' ones(1,Nr)*'C']); %types of variables
MILPproblem.x0=[]; %initial solution
if exist('param','var')
     info= solveCobraMILP(MILPproblem,param);
else
     info= solveCobraMILP(MILPproblem);
end
if SM >= 3
    disp(info.origStat);
end
Cor=[];
if isempty(info.full)
    yn=true;
    if SM >= 3
        disp('Non-decomposable1. EFM');
    end
    p2=[];
    Rz=[];
    Cor=[];
else
    p2=info.cont;
    %p2=sol(Nr+1:2*Nr);
    %a1=sol(1:Nr);
    if sum(p2>epsP)>=sum(p2~=0)/3
        nz=(p2>epsP);
        p2(~nz)=0;
        p2=min(flux(nz)./p2(nz))*p2;
        p2(p2<=min(flux(nz))*(10^(-4)))=0;
        imb=sum(abs(ST*p2));
        if sum(nz)<sum(fm>eps8) && sum(nz)>0 && imb<=eps10
            yn=false;
            if SM >= 3
                disp('Decomposable. Not EFM');
                disp({'imbalance of p2:' imb});
            end
            if imb>epsB
                Cor=p2;
                [p2] = fluxcorrect(ST,p2,lb0,flux,'cobra');
                Cor=[Cor p2];
            end
            Rz=~nz&(fm>eps8);
            indr=1:Nr;
            Rz=indr(Rz);
            Rz=Rz(1);
        else
            yn=true;
            if SM >= 3
                disp('Non-decomposable2. EFM');
            end
            p2=[];
            Rz=[];
            Cor=[];
        end
    else
        if sum(p2>eps9)>=2
            p2(p2<=eps9)=0;
            p2=p2/min(p2(p2~=0));
            a1=logical(p2);
            nz=(p2>0.5);
            if sum(a1 & fm>eps8)==sum(a1)
                imb=sum(abs(ST*p2));
                if SM >= 3
                    disp({'imbalance of p2:' imb});
                end
                if imb>eps9
                    yn=true;
                    if SM >= 3
                        disp('Non-decomposable3. EFM');
                    end
                    p2=[];
                    Rz=[];
                    Cor=[];
                else
                    yn=false;
                    if SM >= 3
                        disp('Decomposable. Not EFM');
                    end
                    p2(p2<=eps10)=0;
                    nz2=(p2>0);
                    p2=min(flux(nz2)./p2(nz2))*p2;
                    imb=sum(abs(ST*p2));
                    if SM >= 3
                        disp({'imbalance of p2:' imb});
                    end
                    if imb>epsB
                        Cor=p2;
                        [p2] = fluxcorrect(ST,p2,lb0,flux,'cobra');
                        Cor=[Cor p2];
                    end
                    Rz=~nz&(fm>eps8);
                    indr=1:Nr;
                    Rz=indr(Rz);
                    Rz=Rz(1);
                end
            else
                yn=true;
                if SM >= 3
                    disp('Non-decomposable4. EFM');
                end
                p2=[];
                Rz=[];
                Cor=[];
            end
        else
            yn=true;
            if SM >= 3
                disp('Non-decomposable4. EFM');
            end
            p2=[];
            Rz=[];
            Cor=[];
        end
    end
end
end

