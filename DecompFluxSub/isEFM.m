function [decomp] = isEFM(FM,ST,lb0,ub0,solver,param)
%[decomp] = isEFM(EFM,ST)
%[decomp] = isEFM(EFM,ST,lb0,ub0)
%Determine whether a given set of flux modes is EFMs or not.
%FM: N by K matrix, K flux modes, each containing N reactions.
%ST: M by N stoichiometrix matrix
%lb0 (optional):lower bound, deafult zeros
%ub0 (optional): upper bound, deafult 1000000
%solver (optional): solver used to solve LP: 'cplex', 'cobra' or 'matlab'
%                   search automatically if not provided.
%param (optional): optimization parameters
%
%decomp: true if column is EFM and false if not.
if nargin==2
    lb0=zeros(size(ST,2),1);
    ub0=1000000*ones(size(ST,2),1);
end
if isempty(lb0)
     lb0=zeros(size(ST,2),1);
end
if isempty(ub0)
    ub0=1000000*ones(size(ST,2),1);
end
if ~exist('solver','var')
    if ~isempty(which('cplexmilp'))
        solver='cplex';
    elseif ~isempty(which('solveCobraMILP'))
        solver='cobra';
    elseif ~isempty(which('linprog'))
        solver='matlab';
    else
        error('A supported solver (Cplex, Cobra or Matlab linprog solver) cannot be identified.')
    end
    disp(['Solver identified: ' solver])
else
    solver=lower(solver);
    if strcmp(solver,'cplex')
        checksolver=which('cplexmilp');
    elseif strcmp(solver,'cobra')
        checksolver=which('solveCobraMILP');
    elseif strcmp(solver,'matlab')
        checksolver=which('linprog');
    else
        checksolver=[];
    end
    if isempty(checksolver)
        error('A supported solver (Cplex or Cobra solver) cannot be identified.')
    end
end

K=size(FM,2);
decomp=false(K,1);
p=find(floor(K/10*(1:10)), 1, 'first');
flux=10000*ones(size(FM,1),1);
for k=1:K
    fm=FM(:,k);
    if exist('param','var')
        [yn,~,~] = isDC(fm,ST,lb0,ub0,flux,k,solver,param);
    else
        [yn,~,~] = isDC(fm,ST,lb0,ub0,flux,k,solver);
    end
    decomp(k)=yn;
    if k==floor(K/10*p)
        while k==floor(K/10*p)
            p=p+1;
        end
        disp([num2str(round(100*k/K)) '% finished.']);
    end
end
end

